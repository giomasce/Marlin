/*
  motion_control.c - high level interface for issuing motion commands
  Part of Grbl

  Copyright (c) 2009-2011 Simen Svale Skogsrud
  Copyright (c) 2011 Sungeun K. Jeon
  Copyright (c) 2016 Giovanni Mascellani <g.mascellani@gmail.com>
  
  Grbl is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Grbl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Grbl.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Marlin.h"
#include "stepper.h"
#include "planner.h"

// The arc is approximated by generating a huge number of tiny, linear segments. The length of each 
// segment is configured in settings.mm_per_arc_segment.  
void mc_arc(float *position, float *target, float *offset, uint8_t axis_0, uint8_t axis_1, 
  uint8_t axis_linear, float feed_rate, float radius, uint8_t isclockwise, uint8_t extruder)
{      
  //   int acceleration_manager_was_enabled = plan_is_acceleration_manager_enabled();
  //   plan_set_acceleration_manager_enabled(false); // disable acceleration management for the duration of the arc
  float center_axis0 = position[axis_0] + offset[axis_0];
  float center_axis1 = position[axis_1] + offset[axis_1];
  float linear_travel = target[axis_linear] - position[axis_linear];
  float extruder_travel = target[E_AXIS] - position[E_AXIS];
  float r_axis0 = -offset[axis_0];  // Radius vector from center to current location
  float r_axis1 = -offset[axis_1];
  float rt_axis0 = target[axis_0] - center_axis0;
  float rt_axis1 = target[axis_1] - center_axis1;
  
  // CCW angle between position and target from circle center. Only one atan2() trig computation required.
  float angular_travel = atan2(r_axis0*rt_axis1-r_axis1*rt_axis0, r_axis0*rt_axis0+r_axis1*rt_axis1);
  if (angular_travel < 0) { angular_travel += 2*M_PI; }
  if (isclockwise) { angular_travel -= 2*M_PI; }
  
  //20141002:full circle for G03 did not work, e.g. G03 X80 Y80 I20 J0 F2000 is giving an Angle of zero so head is not moving
  //to compensate when start pos = target pos && angle is zero -> angle = 2Pi
  if (position[axis_0] == target[axis_0] && position[axis_1] == target[axis_1] && angular_travel == 0)
  {
	  angular_travel += 2*M_PI;
  }
  //end fix G03
  
  float millimeters_of_travel = hypot(angular_travel*radius, fabs(linear_travel));
  if (millimeters_of_travel < 0.001) { return; }
  uint16_t segments = floor(millimeters_of_travel/MM_PER_ARC_SEGMENT);
  if(segments == 0) segments = 1;
  
  /*  
    // Multiply inverse feed_rate to compensate for the fact that this movement is approximated
    // by a number of discrete segments. The inverse feed_rate should be correct for the sum of 
    // all segments.
    if (invert_feed_rate) { feed_rate *= segments; }
  */
  float theta_per_segment = angular_travel/segments;
  float linear_per_segment = linear_travel/segments;
  float extruder_per_segment = extruder_travel/segments;
  
  /* Vector rotation by transformation matrix: r is the original vector, r_T is the rotated vector,
     and phi is the angle of rotation. Based on the solution approach by Jens Geisler.
         r_T = [cos(phi) -sin(phi);
                sin(phi)  cos(phi] * r ;
     
     For arc generation, the center of the circle is the axis of rotation and the radius vector is 
     defined from the circle center to the initial position. Each line segment is formed by successive
     vector rotations. This requires only two cos() and sin() computations to form the rotation
     matrix for the duration of the entire arc. Error may accumulate from numerical round-off, since
     all double numbers are single precision on the Arduino. (True double precision will not have
     round off issues for CNC applications.) Single precision error can accumulate to be greater than
     tool precision in some cases. Therefore, arc path correction is implemented. 

     Small angle approximation may be used to reduce computation overhead further. This approximation
     holds for everything, but very small circles and large mm_per_arc_segment values. In other words,
     theta_per_segment would need to be greater than 0.1 rad and N_ARC_CORRECTION would need to be large
     to cause an appreciable drift error. N_ARC_CORRECTION~=25 is more than small enough to correct for 
     numerical drift error. N_ARC_CORRECTION may be on the order a hundred(s) before error becomes an
     issue for CNC machines with the single precision Arduino calculations.
     
     This approximation also allows mc_arc to immediately insert a line segment into the planner 
     without the initial overhead of computing cos() or sin(). By the time the arc needs to be applied
     a correction, the planner should have caught up to the lag caused by the initial mc_arc overhead. 
     This is important when there are successive arc motions. 
  */
  // Vector rotation matrix values
  float cos_T = 1-0.5*theta_per_segment*theta_per_segment; // Small angle approximation
  float sin_T = theta_per_segment;
  
  float arc_target[4];
  float sin_Ti;
  float cos_Ti;
  float r_axisi;
  uint16_t i;
  int8_t count = 0;

  // Initialize the linear axis
  arc_target[axis_linear] = position[axis_linear];
  
  // Initialize the extruder axis
  arc_target[E_AXIS] = position[E_AXIS];

  for (i = 1; i<segments; i++) { // Increment (segments-1)
    
    if (count < N_ARC_CORRECTION) {
      // Apply vector rotation matrix 
      r_axisi = r_axis0*sin_T + r_axis1*cos_T;
      r_axis0 = r_axis0*cos_T - r_axis1*sin_T;
      r_axis1 = r_axisi;
      count++;
    } else {
      // Arc correction to radius vector. Computed only every N_ARC_CORRECTION increments.
      // Compute exact location by applying transformation matrix from initial radius vector(=-offset).
      cos_Ti = cos(i*theta_per_segment);
      sin_Ti = sin(i*theta_per_segment);
      r_axis0 = -offset[axis_0]*cos_Ti + offset[axis_1]*sin_Ti;
      r_axis1 = -offset[axis_0]*sin_Ti - offset[axis_1]*cos_Ti;
      count = 0;
    }

    // Update arc_target location
    arc_target[axis_0] = center_axis0 + r_axis0;
    arc_target[axis_1] = center_axis1 + r_axis1;
    arc_target[axis_linear] += linear_per_segment;
    arc_target[E_AXIS] += extruder_per_segment;

    clamp_to_software_endstops(arc_target);
    plan_buffer_line(arc_target[X_AXIS], arc_target[Y_AXIS], arc_target[Z_AXIS], arc_target[E_AXIS], feed_rate, extruder);
    
  }
  // Ensure last segment arrives at target location.
  plan_buffer_line(target[X_AXIS], target[Y_AXIS], target[Z_AXIS], target[E_AXIS], feed_rate, extruder);

  //   plan_set_acceleration_manager_enabled(acceleration_manager_was_enabled);
}

// See the meaning in the documentation of mc_cubic().
#define MIN_STEP 0.002
#define MAX_STEP 0.1
#define SIGMA 0.1

/* Compute the linear interpolation between to real numbers.
*/
inline static float interp(float a, float b, float t) {
  return (1.0 - t) * a + t * b;
}

/* Compute a Bézier curve using the De Casteljau's algorithm (see
   https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm), which is
   easy to code and has good numerical stability (very important,
   since Arudino works with limited precision real numbers).
 */
inline static float eval_bezier(float a, float b, float c, float d, float t) {
  float iab = interp(a, b, t);
  float ibc = interp(b, c, t);
  float icd = interp(c, d, t);
  float iabc = interp(iab, ibc, t);
  float ibcd = interp(ibc, icd, t);
  float iabcd = interp(iabc, ibcd, t);
  return iabcd;
}

/* We approximate Euclidean distance with the sum of the coordinates
   offset (so-called "norm 1"), which is quicker to compute.
 */
inline static float dist1(float x1, float y1, float x2, float y2) {
  return fabs(x1 - x2) + fabs(y1 - y2);
}

/* The algorithm for computing the step is loosely based on the one in
   Kig (see
   https://sources.debian.net/src/kig/4:15.08.3-1/misc/kigpainter.cpp/#L759);
   however, we do not use the stack.

   The algorithm goes as it follows: the parameters t runs from 0.0 to
   1.0 describing the curve, which is evaluated by eval_bezier(). At
   each iteration we have to choose a step, i.e., the increment of the
   t variable. By default the step of the previous iteration is taken,
   and then it is enlarged or reduced depending on how straight the
   curve locally is. The step is always clamped between MIN_STEP/2 and
   2*MAX_STEP. MAX_STEP is taken at the first iteration.

   For some t, the step value is considered acceptable if the curve in
   the interval [t, t+step] is sufficiently straight, i.e.,
   sufficiently close to linear interpolation. In practice the
   following test is performed: the distance between eval_bezier(...,
   t+step/2) is evaluated and compared with 0.5*(eval_bezier(...,
   t)+eval_bezier(..., t+step)). If it is smaller than SIGMA, then the
   step value is considered acceptable, otherwise it is not. The code
   seeks to find the larger step value which is considered acceptable.

   At every iteration the recorded step value is considered and then
   iteratively halved until it becomes acceptable. If it was already
   acceptable in the beginning (i.e., no halving were done), then
   maybe it was necessary to enlarge it; then it is iteratively
   doubled while it remains acceptable. The last acceptable value
   found is taken, provided that it is between MIN_STEP and MAX_STEP
   and does not bring t over 1.0.

   Caveat: this algorithm is not perfect, since it can happen that a
   step is considered acceptable even when the curve is not linear at
   all in the interval [t, t+step] (but its mid point coincides "by
   chance" with the midpoint according to the parametrization). This
   kind of glitches can be eliminated with proper first derivative
   estimates; however, given the improbability of such configurations,
   the mitigation offered by MIN_STEP and the small computational
   power available on Arduino, I think it is not wise to implement it.
 */
void mc_cubic(float *position, float *target, float *offset, uint8_t axis_0, uint8_t axis_1, uint8_t axis_linear, float feed_rate, uint8_t extruder) {
  // Absolute first and second control points are recovered.
  float first0 = position[axis_0] + offset[0];
  float first1 = position[axis_1] + offset[1];
  float second0 = target[axis_0] + offset[2];
  float second1 = target[axis_1] + offset[3];
  float t = 0.0;

  float tmp[4];
  tmp[axis_0] = position[axis_0];
  tmp[axis_1] = position[axis_1];
  float step = MAX_STEP;
  while (t < 1.0) {
    // First try to reduce the step in order to make it sufficiently
    // close to a linear interpolation.
    bool did_reduce = false;
    float new_t = t + step;
    if (new_t > 1.0) {
      new_t = 1.0;
    }
    float new_pos0 = eval_bezier(position[axis_0], first0, second0, target[axis_0], new_t);
    float new_pos1 = eval_bezier(position[axis_1], first1, second1, target[axis_1], new_t);
    while (true) {
      if (new_t - t < MIN_STEP) {
        break;
      }
      float candidate_t = 0.5 * (t + new_t);
      float candidate_pos0 = eval_bezier(position[axis_0], first0, second0, target[axis_0], candidate_t);
      float candidate_pos1 = eval_bezier(position[axis_1], first1, second1, target[axis_1], candidate_t);
      float interp_pos0 = 0.5 * (tmp[axis_0] + new_pos0);
      float interp_pos1 = 0.5 * (tmp[axis_1] + new_pos1);
      if (dist1(candidate_pos0, candidate_pos1, interp_pos0, interp_pos1) > SIGMA) {
        new_t = candidate_t;
        new_pos0 = candidate_pos0;
        new_pos1 = candidate_pos1;
        did_reduce = true;
      } else {
        break;
      }
    }

    // If we did not reduce the step, maybe we should enlarge it.
    if (!did_reduce) {
      while (true) {
        if (new_t - t > MAX_STEP) {
          break;
        }
        float candidate_t = t + 2.0 * (new_t - t);
        if (candidate_t > 1.0) {
          break;
        }
        float candidate_pos0 = eval_bezier(position[axis_0], first0, second0, target[axis_0], candidate_t);
        float candidate_pos1 = eval_bezier(position[axis_1], first1, second1, target[axis_1], candidate_t);
        float interp_pos0 = 0.5 * (tmp[axis_0] + candidate_pos0);
        float interp_pos1 = 0.5 * (tmp[axis_1] + candidate_pos1);
        if (dist1(new_pos0, new_pos1, interp_pos0, interp_pos1) > SIGMA) {
          break;
        } else {
          new_t = candidate_t;
          new_pos0 = candidate_pos0;
          new_pos1 = candidate_pos1;
        }
      }
    }

    // Check some postcondition; they are disabled in the actual
    // Marlin build, but if you test the same code on a computer you
    // may want to check they are respect.
    //assert(new_t <= 1.0);
    //assert(new_t - t >= MIN_STEP / 2.0);
    //assert(new_t - t <= MAX_STEP * 2.0);

    step = new_t - t;
    t = new_t;

    // Compute and send new position
    tmp[axis_0] = new_pos0;
    tmp[axis_1] = new_pos1;
    // FIXME. The following two are wrong, since the parameter t is
    // not linear in the distance.
    tmp[axis_linear] = interp(position[axis_linear], target[axis_linear], t);
    tmp[E_AXIS] = interp(position[E_AXIS], target[E_AXIS], t);
    clamp_to_software_endstops(tmp);
    plan_buffer_line(tmp[X_AXIS], tmp[Y_AXIS], tmp[Z_AXIS], tmp[E_AXIS], feed_rate, extruder);
  }
}
