
#include <t8_cmesh.h>
#include <t8_vec.h>
#include <t8_cmesh_vtk.h>
#include <t8dg.h>
#include <t8dg_mptrac.h>
#include "t8dg_common.h"
#include "../t8/example/mptrac/t8_mptrac_interpolate.h"
#include <math.h>

t8dg_scalar_function_3d_time_fn
t8dg_common_initial_cond_fn (int initial_cond_arg)
{
  switch (initial_cond_arg) {
  case (0):
    return t8dg_scalar3d_constant_one;
  case (1):
    return t8dg_scalar1d_hat_function;
  case (2):
    return t8dg_scalar1d_step_function;
  case (3):
    return t8dg_scalar3d_cos_product;
  case (4):
    return t8dg_scalar3d_norm_function;
  case (5):
    return t8dg_scalar2d_hat_function;
  case (6):
    return t8dg_scalar2d_step_function;
  case (7):
    return t8dg_scalar2d_triangle_step_function;
  case (8):
    return t8dg_scalar3d_step_function;
  case (9):
    return t8dg_circle_ring_step_function;
  case (10):
    return t8dg_scalar2d_angle;
  case (11):
    return t8dg_cylinder_ring_sin_product_fn;
  case (12):
    return t8dg_cylinder_ring_step_function;
  case (13):
    return t8dg_smooth_indicator1Dfn;
  case (14):
    return t8dg_smooth_indicator2Dfn;
  case (15):
    return t8dg_smooth_indicator3Dfn;
  case (16):
    return t8dg_scalar3d_constant_zero;
  case (17):
    return t8dg_circle_ring_sin_product_fn;
  case (18):
    return t8dg_mptrac_box_indicator_fn;
  case (19):
    return t8dg_williamson_etal_cosine_bell_fn;
  case (20):
    return t8dg_smooth_indicator3D_4Spheres_fn;
  case (21):
    return t8dg_smooth_indicator3D_3Spheres_above_below_fn;
  case (22):
    return t8dg_smooth_indicator3D_bottom_fn;
  case (23):
    return t8dg_smooth_indicator3D_4Spheres_between_fn;
  case (24):
    return t8dg_smooth_indicator3D_boundary_x_fn;
  default:
    return NULL;
  }
}

t8dg_scalar_function_3d_time_fn
t8dg_common_analytic_solution_fn (int initial_cond_arg, double diffusion_coefficient)
{
  if (diffusion_coefficient > 0) {
    switch (initial_cond_arg) {
    case (0):
      return t8dg_scalar3d_constant_one;
    case (3):
      return t8dg_scalar3d_cos_product;
    default:
      return NULL;
    }
  }
  else {
    return t8dg_common_initial_cond_fn (initial_cond_arg);
  }
}

double
t8dg_scalar3d_constant_one (const double x[3], const double t, void *fn_data)
{
  return 1;
}

double
t8dg_scalar2d_angle (const double x[3], const double t, void *fn_data)
{
  return fabs (atan2 (x[1], x[0]));
}

double
t8dg_scalar3d_cos_product (const double x[3], const double t, void *fn_data)
{
  t8dg_scalar3d_cos_product_data_t *cos_data = (t8dg_scalar3d_cos_product_data_t *) fn_data;
  int                 dimension = cos_data->dim;
  double              diffusion_coefficient = cos_data->diffusion_coefficient;
  return exp (-diffusion_coefficient * dimension * 4 * M_PI * M_PI * t) * cos (2 * M_PI * x[0]) * cos (2 * M_PI * x[1]) * cos (2 * M_PI *
                                                                                                                               x[2]);
}

double
t8dg_scalar3d_constant_zero (const double x[3], const double t, void *fn_data)
{
  return 0;
}

double
t8dg_scalar1d_hat_function (const double x[3], const double t, void *fn_data)
{
  return 0.5 - (fabs (0.5 - x[0]));
}

double
t8dg_scalar2d_hat_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return sqrt (0.5) - t8_vec_dist (x, center);
}

double
t8dg_scalar3d_norm_function (const double x[3], const double t, void *fn_data)
{
  return t8_vec_norm (x);
}

double
t8dg_scalar1d_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0, 0 };
  return t8_vec_dist (x, center) < 0.2;
}

double
t8dg_scalar2d_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0.5, 0 };
  return t8_vec_dist (x, center) < 0.15;
}

double
t8dg_scalar2d_triangle_step_function (const double x[3], const double t, void *fn_data)
{
  return x[0] > 0.3 && x[1] > 0.3 && x[0] + x[1] < 0.9;
}

double
t8dg_scalar3d_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 0.5, 0.5, 0.5 };
  return t8_vec_dist (x, center) < 0.15;
}

double
t8dg_circle_ring_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 1.5, 0.0, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < 0.1)
    return 1;
  if (dist > 0.2)
    return 0;
  dist = (dist - 0.1) * 10;     /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_cylinder_ring_sin_product_fn (const double x[3], const double t, void *fn_data)
{
  double              angle = atan2 (x[1], x[0]);
  double              radius = sqrt (x[0] * x[0] + x[1] * x[1]);
  double              h = x[2];
  return sin (angle) * sin ((radius - 1.5) * 2 * M_PI) * sin ((h - 0.5) * 2 * M_PI);
}

double
t8dg_circle_ring_sin_product_fn (const double x[3], const double t, void *fn_data)
{
  double              angle = atan2 (x[1], x[0]);
  double              radius = sqrt (x[0] * x[0] + x[1] * x[1]);
  return sin (angle) * sin ((radius - 1.5) * 2 * M_PI);
}

static double
t8dg_smooth_h (const double x)
{
  if (x <= 0)
    return 0;
  return exp (-1 / x);
}

double
t8dg_smooth_g (const double x)
{
  return t8dg_smooth_h (1 - x) / (t8dg_smooth_h (x) + t8dg_smooth_h (1 - x));
}

double
t8dg_cylinder_ring_step_function (const double x[3], const double t, void *fn_data)
{
  double              center[3] = { 1.5, 0.0, 0.5 };
  double              dist = t8_vec_dist (x, center);
  if (dist < 0.1)
    return 1;
  if (dist > 0.2)
    return 0;
  dist = (dist - 0.1) * 10;     /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_cylinder_ring_source_fn (const double x[3], const double t, void *fn_data)
{
  return 30 * t8dg_cylinder_ring_step_function (x, t, fn_data);
}

double
t8dg_cos_indicator1Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.1;
  double              center[3] = { 0.5, 0.0, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_cos_indicator2Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.1;
  double              center[3] = { 0.5, 0.5, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_cos_indicator3Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.1;
  double              center[3] = { 0.5, 0.5, 0.5 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return (cos (dist * M_PI) + 1) / 2;
}

double
t8dg_smooth_indicator1Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.5;
  double              center[3] = { 0.5, 0.0, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_smooth_indicator2Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 0.2;
  double              smoothing_factor = 0.5;
  double              center[3] = { 0.5, 0.5, 0.0 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_smooth_indicator3Dfn (const double x[3], const double t, void *fn_data)
{
  double              radius = 1. / (2 * 360);  // 0.5 deg
  double              smoothing_factor = 50;    // 2.5 deg
  double              center[3] = { 0.5, 0.5, 0.5 };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_smooth_indicator3D_bottom_fn (const double x[3], const double t, void *fn_data)
{
  double              radius = 1. / (2 * 360);
  double              smoothing_factor = 40;
  double              speed_in_kmph = 1.0;      //TODO: Read out of windfile
  double              max_altitude_in_km = 60.0;        //TODO: Read out of windfile
  double              center_z = 0.1 + (speed_in_kmph * t) / max_altitude_in_km;

  double              center[3] = { 0.5, 0.5, center_z };
  double              dist = t8_vec_dist (x, center);
  if (dist < radius)
    return 1;
  if (dist > (1 + smoothing_factor) * radius)
    return 0;
  dist = (dist - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist);
}

double
t8dg_smooth_indicator3D_boundary_x_fn (const double x[3], const double t, void *fn_data)
{
  double              radius = 1. / (2 * 360);
  double              smoothing_factor = 40;

  double              center[3] = { 0, 0.5, 0.5 };
  /* Array for x coordinates of different center for spheres */
  /* Center x=0 and x=1 coincide because of periodic boundaries */
  double              center_x[2] = { 0., 1. };
  int                 center_x_len = sizeof (center_x) / sizeof (double);

  double              dist;
  double              dist_min = t8_vec_dist (x, center);

  for (int idx = 0; idx < center_x_len; idx++) {
    center[0] = center_x[idx];
    dist = t8_vec_dist (x, center);

    if (dist < dist_min)
      dist_min = dist;

    if (dist < radius)
      return 1;
  }

  if (dist_min > (1 + smoothing_factor) * radius)
    return 0;

  dist_min = (dist_min - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist_min);
}

double
t8dg_smooth_indicator3D_4Spheres_fn (const double x[3], const double t, void *fn_data)
{
  double              radius = 1. / (2 * 360);  // 0.5 deg
  double              smoothing_factor = 50;    // 2.5 deg

  double              center[3] = { 0, 0.5, 0.5 };
  /* Array for x coordinates of different center for spheres */
  /* Center x=0 and x=1 coincide because of periodic boundaries */
  double              center_x[5] = { 0, 0.25, 0.5, 0.75, 1 };
  int                 center_x_len = sizeof (center_x) / sizeof (double);

  double              dist;
  double              dist_min = t8_vec_dist (x, center);

  for (int idx = 0; idx < center_x_len; idx++) {
    center[0] = center_x[idx];
    dist = t8_vec_dist (x, center);

    if (dist < dist_min)
      dist_min = dist;

    if (dist < radius)
      return 1;
  }

  if (dist_min > (1 + smoothing_factor) * radius)
    return 0;

  dist_min = (dist_min - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist_min);
}

double
t8dg_smooth_indicator3D_4Spheres_between_fn (const double x[3], const double t, void *fn_data)
{
  double              radius = 1. / (2 * 360);  // 0.5 deg
  double              smoothing_factor = 50;    // 2.5 deg

  double              center[3] = { 0.125, 0.5, 0.5 };
  /* Array for x coordinates of different center for spheres */
  /* Center x=0 and x=1 coincide because of periodic boundaries */
  double              center_x[4] = { 0.125, 0.375, 0.625, 0.875 };
  int                 center_x_len = sizeof (center_x) / sizeof (double);

  double              dist;
  double              dist_min = t8_vec_dist (x, center);

  for (int idx = 0; idx < center_x_len; idx++) {
    center[0] = center_x[idx];
    dist = t8_vec_dist (x, center);

    if (dist < dist_min)
      dist_min = dist;

    if (dist < radius)
      return 1;
  }

  if (dist_min > (1 + smoothing_factor) * radius)
    return 0;

  dist_min = (dist_min - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist_min);
}

double
t8dg_smooth_indicator3D_3Spheres_above_below_fn (const double x[3], const double t, void *fn_data)
{
  double              radius = 1. / (2 * 360);  // 0.5 deg
  double              smoothing_factor = 50;    // 2.5 deg

  double              center[3] = { 0.5, 0.25, 0.5 };
  /* Array for x coordinates of different center for spheres */
  double              center_y[3] = { 0.25, 0.5, 0.75 };
  int                 center_y_len = sizeof (center_y) / sizeof (double);

  double              dist;
  double              dist_min = t8_vec_dist (x, center);

  for (int idx = 0; idx < center_y_len; idx++) {
    center[1] = center_y[idx];
    dist = t8_vec_dist (x, center);

    if (dist < dist_min)
      dist_min = dist;

    if (dist < radius)
      return 1;
  }

  if (dist_min > (1 + smoothing_factor) * radius)
    return 0;

  dist_min = (dist_min - radius) / (radius * smoothing_factor); /* transform to [0,1] */
  return t8dg_smooth_g (dist_min);
}

double
t8dg_williamson_etal_cosine_bell_fn (const double x[3], const double t, void *fn_data)
{
  /* Transformed test initial condition by Williamson et al., 1992; Test Case 1, p. 217 */
  /* Original version was made for 2D SWE  */
  /* Idea to transform to 3D:
   *    - Transform z coordinate of cube \in [0,1] to altitude via interpolation
   *    - To do so, transform min and max pressure to max and min altitude
   *    - Build up an 'ellipse' in z direction shrank by a (arbitrary) factor
   */

  /* Calculations are done in km as well as all length values */

  t8dg_williamson_etal_data_t *williamson_data = (t8dg_williamson_etal_data_t *) fn_data;
  double              inner_radius_factor = williamson_data->inner_radius_factor;
  double              ring_radius_factor = williamson_data->ring_radius_factor;
  double              weight_z_direction_factor = williamson_data->weight_z_direction_factor;

  /* Earth constants  */
  const double        a = 6371.220;

  double              inner_radius = inner_radius_factor * a;
  double              ring_radius = ring_radius_factor * a;
  const double        weight_z_direction = weight_z_direction_factor * a;

  /* Initial values  */
  const int           h_0 = 1;
  double              x_globe[3] = { 0, 0, 0 }; /* needs to be set, transform input x */

  /* Transform deg to rad within interpolation */
  x_globe[0] = ((x[0] - 0.5) * 360.) * M_PI / 180.;
  x_globe[1] = ((x[1] - 0.5) * 180.) * M_PI / 180.;

  /* Found in (arbitrary) output file on JUWELS */
  /* For calculation/interpolation of altitude max and min value Pressure levels: 1013.25, 878.364 ... 0.191952 hPa */
  double              pressure_max = 1013.25;
  double              pressure_min = 0.191952;
  double              altitude_min_in_km = Z (pressure_max);
  double              altitude_max_in_km = Z (pressure_min);
  double              altitude_middle_in_km = (altitude_max_in_km + altitude_min_in_km) * 0.5;

  /* The paper gives (3/2*pi, 0) as 2D center. We set it to (0.5,0.5,0.5) on cube. */
  /* Lon and lat needs to be set to 0 to place center in the middle of cube. */
  const double        center[3] = { 0., 0., altitude_middle_in_km };

  /* Great circle distance between input (transformed to globe) and center */
  double              r = a * acos (sin (center[1]) * sin (x_globe[1]) + cos (center[1]) * cos (x_globe[1]) * cos (center[0] - x_globe[0]));

  /* Convert altitude in cube to x_globe[2] [km] via scaling. */
  x_globe[2] = (1 - x[2]) * altitude_min_in_km + x[2] * altitude_max_in_km;

  double              z_dist = fabs (x_globe[2] - center[2]);

  /* Value of initial function; depending on x_globe and its great circle distance */
  double              concentration = 0;

  /* Ellipse around center to control expansion in z direction; dist_z=radius overshoots the cube and creates unphysical initial condition */
  /* z_dist is weighted with sqrt (weight_z_direction) */
  double              dist = sqrt (r * r + weight_z_direction * z_dist * z_dist);

  if (dist < inner_radius)
    return 1;
  if (dist > ring_radius + inner_radius)
    return 0;

  dist = (dist - inner_radius) / ring_radius; /* transform to [0,1] */

  return (h_0 * 0.5) * (1 + cos (M_PI * dist));
}
