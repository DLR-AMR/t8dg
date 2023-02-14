#include <t8_forest.h>
#include "t8dg_dof.h"
#include "t8dg_values.h"
#include "t8dg_adapt.h"
#include <t8_element_cxx.hxx>
#include <t8_vec.h>
#include <t8dg_mptrac.h>

t8_forest_adapt_t
t8dg_adapt_fn_arg (int adapt_arg)
{
  switch (adapt_arg) {
  case 0:
    return t8dg_adapt_mass;
    break;
  case 1:
    return t8dg_adapt_rel_min_max;
    break;
  case 2:
    return t8dg_adapt_smooth_indicator;
    break;
  case 3:
    return t8dg_adapt_smooth_indicator_hypercube;
    break;
  case 4:
    return t8dg_adapt_mptrac_hypercube;
    break;

  default:
    T8DG_ABORT ("Wrong adapt fn arg");
    break;
  }
}

t8dg_adapt_data_t  *
t8dg_adapt_data_new (t8dg_values_t * dg_values, int initial_level, int min_level, int max_level, int adapt_fn_arg,
                     int adapt_freq, t8dg_scalar_function_3d_time_fn source_sink_fn, void *source_sink_data,
                     double refinement_threshold, double coarsening_threshold)
{
  t8dg_adapt_data_t  *adapt_data;
  adapt_data = T8DG_ALLOC_ZERO (t8dg_adapt_data_t, 1);
  adapt_data->dg_values = dg_values;
  adapt_data->minimum_refinement_level = min_level;
  adapt_data->initial_refinement_level = initial_level;
  adapt_data->maximum_refinement_level = max_level;
  adapt_data->adapt_fn_arg = adapt_fn_arg;
  adapt_data->adapt_fn = t8dg_adapt_fn_arg (adapt_fn_arg);
  adapt_data->adapt_freq = adapt_freq;
  if (min_level == max_level) {
    adapt_data->adapt_freq = 0;
  }
  adapt_data->source_sink_fn = source_sink_fn;
  adapt_data->source_sink_data = source_sink_data;
  adapt_data->dim = t8dg_values_get_dim (dg_values);
  adapt_data->refinement_threshold = refinement_threshold;
  adapt_data->coarsening_threshold = coarsening_threshold;
  return adapt_data;
}

void
t8dg_adapt_data_interpolate_source_fn (t8dg_adapt_data_t * adapt_data)
{
  double              time = 0;
   /*TODO*/
    adapt_data->source_sink_dof =
    t8dg_dof_values_new (t8dg_values_get_forest (adapt_data->dg_values), t8dg_values_get_global_values_array (adapt_data->dg_values));
  t8dg_values_interpolate_scalar_function_3d_time (adapt_data->dg_values, adapt_data->source_sink_fn, time, adapt_data->source_sink_data,
                                                   adapt_data->source_sink_dof);
}

void
t8dg_adapt_data_destroy (t8dg_adapt_data_t ** p_adapt_data)
{
  t8dg_adapt_data_t  *adapt_data;
  adapt_data = *p_adapt_data;
  T8DG_FREE (adapt_data);
  p_adapt_data = NULL;
}

void
t8dg_adapt_data_set_time (t8dg_adapt_data_t * adapt_data, double time)
{
  adapt_data->time = time;
}

int
t8dg_adapt_uniform (t8_forest_t forest,
                    t8_forest_t forest_from,
                    t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  int                 level;
  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);
  level = ts->t8_element_level (elements[0]);
  return level < adapt_data->initial_refinement_level;
}

int
t8dg_adapt_mass (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  sc_array_t         *element_dof;
  int                 level;
  double              norm, area;
  int                 ifamilyelement;
  double              refinement_threshold = 0.01;
  double              coarsening_threshold = 0.0025;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);
  level = ts->t8_element_level (elements[0]);

  if (level == adapt_data->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }

  element_dof = t8dg_dof_values_new_element_dof_values_view (adapt_data->dof_values, itree, ielement);

  norm = t8dg_values_element_norm_l2_squared (adapt_data->dg_values, element_dof, itree, ielement);
  area = t8dg_values_element_area (adapt_data->dg_values, itree, ielement);

  sc_array_destroy (element_dof);

  if (level < adapt_data->maximum_refinement_level) {

    if (norm / area > refinement_threshold) {
      return 1;
    }
    if (adapt_data->source_sink_fn != NULL) {
      //Check source and sink function
      double              max_value, min_value;
      max_value = t8dg_dof_values_get_max_value (adapt_data->source_sink_dof, itree, ielement);
      min_value = t8dg_dof_values_get_min_value (adapt_data->source_sink_dof, itree, ielement);
      if (min_value != 0 || max_value != 0) {
        return 1;
      }
    }

  }

  if (num_elements > 1) {
    if (level == adapt_data->minimum_refinement_level) {
      return 0;                 /* It is not possible to coarsen this element. If this is wanted, balance is needed outside */
    }

    for (ifamilyelement = 0; ifamilyelement < num_elements; ifamilyelement++) {
      element_dof = t8dg_dof_values_new_element_dof_values_view (adapt_data->dof_values, itree, ielement + ifamilyelement);
      norm = t8dg_values_element_norm_l2_squared (adapt_data->dg_values, element_dof, itree, ielement + ifamilyelement);
      area = t8dg_values_element_area (adapt_data->dg_values, itree, ielement + ifamilyelement);
      sc_array_destroy (element_dof);

      if (norm / area > coarsening_threshold) {
        return 0;
      }
    }
    return -1;
  }
  return 0;
}

int
t8dg_adapt_smooth_indicator (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  sc_array_t         *element_dof;
  int                 level;
  double              norm, area, ratio;
  int                 ifamilyelement;
  double              threshold = 0.00000005;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);
  level = ts->t8_element_level (elements[0]);

  if (level == adapt_data->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }

  element_dof = t8dg_dof_values_new_element_dof_values_view (adapt_data->dof_values, itree, ielement);

  norm = t8dg_values_element_norm_l2_squared (adapt_data->dg_values, element_dof, itree, ielement);
  area = t8dg_values_element_area (adapt_data->dg_values, itree, ielement);
  ratio = norm / area;

  sc_array_destroy (element_dof);

  if (ratio > 2 * threshold && ratio < 1 - 2 * threshold) {
    return level < adapt_data->maximum_refinement_level;
  }

  if (num_elements > 1) {
    if (level == adapt_data->minimum_refinement_level) {
      return 0;                 /* It is not possible to coarsen this element. If this is wanted, balance is needed outside */
    }

    for (ifamilyelement = 0; ifamilyelement < num_elements; ifamilyelement++) {
      element_dof = t8dg_dof_values_new_element_dof_values_view (adapt_data->dof_values, itree, ielement + ifamilyelement);
      norm = t8dg_values_element_norm_l2_squared (adapt_data->dg_values, element_dof, itree, ielement + ifamilyelement);
      area = t8dg_values_element_area (adapt_data->dg_values, itree, ielement + ifamilyelement);
      ratio = norm / area;
      sc_array_destroy (element_dof);

      if (ratio > threshold && ratio < 1 - threshold) {
        return 0;
      }
    }
    return -1;
  }
  return 0;
}

int
t8dg_adapt_mptrac_hypercube (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements,
                             t8_element_t * elements[])
{
  const t8_element_t *element = elements[0];
  const t8dg_adapt_data_t  *adapt_data = (const t8dg_adapt_data_t *) t8_forest_get_user_data (forest);
  const int level = ts->t8_element_level (elements[0]);

#if 0
  double midpoint[3];
  double corner_0[3];
  double corner_7[3];

  /* Check whether corner 0, midpoint or corner 7 is contained in the box source.
   * If so, refine the element. */
  t8_forest_element_centroid (forest_from, itree, element, midpoint);
  t8_forest_element_coordinate (forest_from, itree, element, 0, corner_0);
  t8_forest_element_coordinate (forest_from, itree, element, 7, corner_7);
  if (t8dg_mptrac_box_source (midpoint, -1, adapt_data->source_sink_data)
  || t8dg_mptrac_box_source (corner_0, -1, adapt_data->source_sink_data)
  || t8dg_mptrac_box_source (corner_7, -1, adapt_data->source_sink_data)) {
    return 1 && level < adapt_data->maximum_refinement_level;
  }
  /* If the above does not trigger, the box source may be too small to be detected
   * by an element. Therefore, we also refine if the lower left corner of the box lies inside
   * this element. */
  const double box_point[3] = {0.5, 0.5, 0.5};
  if (t8_forest_element_point_inside (forest_from, itree, element, box_point, 1e-3)) {
    return 1 && level < adapt_data->maximum_refinement_level;
  }
#endif

#if 0
  /* Refine at top/bottom and poles */

  int is_upper_boundary = 0;
  int is_lower_boundary = 0;
  int is_south_boundary = 0;
  int is_north_boundary = 0;
  /* Check for upper boundary of unit cube */
  if (ts->t8_element_is_root_boundary (element, 5)) {
    is_upper_boundary = 1;
  }
  /* Check for upper boundary of unit cube */
  if (ts->t8_element_is_root_boundary (element, 4)) {
    is_lower_boundary = 1;
  }
  /* Check for south boundary of unit cube */
  if (ts->t8_element_is_root_boundary (element, 2)) {
    is_south_boundary = 1;
  }
  /* Check for north boundary of unit cube */
  if (ts->t8_element_is_root_boundary (element, 3)) {
    is_north_boundary = 1;
  }
  const int is_at_pole = is_south_boundary || is_north_boundary;
  const int is_at_top_or_bottom = is_upper_boundary || is_lower_boundary;


  if ((is_at_top_or_bottom || is_at_pole) && level < adapt_data->maximum_refinement_level) {
    return 1;
  }
#endif

  /* In the very first time step we refine around our source box. */
  if (adapt_data->time == 0) {
    const t8dg_mptrac_flux_data* flux_data = (const t8dg_mptrac_flux_data*) adapt_data->source_sink_data;

    double midpoint[3];
    double vert_dist_km, horiz_dist_deg;
    const double vertical_distance_threshold_km = 400; /* Is halfed with each refinement level */
    const double horiz_distance_threshold_deg = 1200; /* Is halfed with each refinement level */
    const double max_vert_distance = vertical_distance_threshold_km / (1 << level);
    const double max_horiz_distance = horiz_distance_threshold_deg / (1 << level);

    /* Compute midpoint and level of the current element. */
    t8_forest_element_centroid (forest_from, itree, element, midpoint);

    /* Determine the distance of the midpoint from the box. */  
    flux_data->point_max_vert_and_horiz_distance_from_box (midpoint, &vert_dist_km, &horiz_dist_deg);

    t8dg_debugf ("level %i distance_vert %f distance_horiz %f max_dist (%f, %f)\n", level, vert_dist_km, horiz_dist_deg, max_vert_distance, max_horiz_distance);
    if (vert_dist_km < max_vert_distance && horiz_dist_deg < max_horiz_distance) {
      return 1 && level < adapt_data->maximum_refinement_level;
    }
  }

  return t8dg_adapt_smooth_indicator (forest, forest_from, itree, ielement, ts, num_elements, elements);
#if 0
  /* Pass on to the rel min/max criterion to refine according to gradient. */
  return t8dg_adapt_rel_min_max (forest, forest_from, itree, ielement, ts, num_elements, elements);
#endif
}                                       

int
t8dg_adapt_smooth_indicator_hypercube (t8_forest_t forest,
                                       t8_forest_t forest_from,
                                       t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements,
                                       t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  int                 level;
  double              element_midpoint[3];
  int                 ifamilyelement;

  double              radius = 0.2, smoothing_factor = 0.5;
  double              indicator_midpoint[3] = { 0, 0, 0 };
  double              indicator_left_midpoint[3] = { 0, 0, 0 };
  double              indicator_right_midpoint[3] = { 0, 0, 0 };
  double              dist, dist_left, dist_middle, dist_right;
  double              threshold = 0.01;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);
  level = ts->t8_element_level (elements[0]);

  if (level == adapt_data->maximum_refinement_level && num_elements == 1) {
    /* It is not possible to refine this level */
    return 0;
  }
  int                 dim = adapt_data->dim;
  int                 idim;
  for (idim = 0; idim < dim; idim++) {
    indicator_midpoint[idim] = 0.5;
    indicator_left_midpoint[idim] = 0.5;
    indicator_right_midpoint[idim] = 0.5;
  }
  indicator_midpoint[0] = 0.5 + adapt_data->time;
  indicator_left_midpoint[0] = indicator_midpoint[0] - 1;
  indicator_right_midpoint[0] = indicator_midpoint[0] + 1;

  t8_forest_element_centroid (forest, itree, elements[0], element_midpoint);

  dist_middle = t8_vec_dist (indicator_midpoint, element_midpoint);
  dist_left = t8_vec_dist (indicator_left_midpoint, element_midpoint);
  dist_right = t8_vec_dist (indicator_right_midpoint, element_midpoint);
  dist = SC_MIN (dist_left, dist_middle);
  dist = SC_MIN (dist, dist_right);

  if (dist > radius - threshold && dist < radius * (1 + smoothing_factor) + threshold) {
    return level < adapt_data->maximum_refinement_level;
  }

  if (num_elements > 1) {
    if (level == adapt_data->minimum_refinement_level) {
      return 0;                 /* It is not possible to coarsen this element. If this is wanted, balance is needed outside */
    }

    for (ifamilyelement = 0; ifamilyelement < num_elements; ifamilyelement++) {
      t8_forest_element_centroid (forest, itree, elements[ifamilyelement], element_midpoint);
      dist = t8_vec_dist (indicator_midpoint, element_midpoint);
      if (dist > radius - threshold && dist < radius * (1 + smoothing_factor) + threshold) {
        return 0;
      }
    }
    return -1;
  }
  return 0;
}

int
t8dg_adapt_rel_min_max (t8_forest_t forest,
                        t8_forest_t forest_from,
                        t8_locidx_t itree, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  int                 level;

  double              max_value;
  double              min_value;
  double              rel_gradient;
  /*Move thresholds into adapt_data? */
  double              gradient_threshold_refine = 0.1;
  double              gradient_threshold_coarsen = 0.01;
  double              min_value_threshold = 0.05;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);

  level = ts->t8_element_level (elements[0]);

  if (level < adapt_data->maximum_refinement_level) {
    max_value = t8dg_dof_values_get_max_value (adapt_data->dof_values, itree, lelement_id);
    min_value = t8dg_dof_values_get_min_value (adapt_data->dof_values, itree, lelement_id);
    if (min_value <= min_value_threshold) {
      if (max_value > min_value_threshold)
        return 1;
    }
    else {
      rel_gradient = (max_value - min_value) / min_value;
      if (rel_gradient > gradient_threshold_refine) {
        return 1;
      }
    }
    if (adapt_data->source_sink_fn != NULL) {
      //Check source and sink function
      max_value = t8dg_dof_values_get_max_value (adapt_data->source_sink_dof, itree, lelement_id);
      min_value = t8dg_dof_values_get_min_value (adapt_data->source_sink_dof, itree, lelement_id);
      if (min_value != 0 || max_value != 0) {
        return 1;
      }
    }
  }

  //not refine, check coarsen
  if (num_elements == 1 || level <= adapt_data->minimum_refinement_level) {
    return 0;
  }
  else {
    int                 ichild;
    for (ichild = 1; ichild < num_elements; ichild++) {
      max_value = t8dg_dof_values_get_max_value (adapt_data->dof_values, itree, lelement_id + ichild);
      min_value = t8dg_dof_values_get_min_value (adapt_data->dof_values, itree, lelement_id + ichild);
      if (min_value <= min_value_threshold && max_value > min_value_threshold)
        return 0;
      rel_gradient = (max_value - min_value) / min_value;
      if (rel_gradient > gradient_threshold_coarsen)
        return 0;
    }
    return -1;
  }
}

void
t8dg_adapt_replace (t8_forest_t forest_old,
                    t8_forest_t forest_new,
                    t8_locidx_t itree,
                    t8_eclass_scheme_c * ts, int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new)
{
  t8dg_adapt_data_t  *adapt_data;
  t8_locidx_t         first_idata_old, first_idata_new;
  int                 ichild, num_children;

  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest_new);
  first_idata_old = t8dg_itree_ielement_to_idata (forest_old, itree, first_ielem_old);
  first_idata_new = t8dg_itree_ielement_to_idata (forest_new, itree, first_ielem_new);
  if (num_elems_old == num_elems_new && num_elems_old == 1) {
    t8dg_values_copy_element_values (adapt_data->dg_values, first_idata_old, first_idata_new);
    t8dg_dof_values_copy_from_index_to_index (adapt_data->dof_values, first_idata_old, adapt_data->dof_values_adapt, first_idata_new);
  }
  else if (num_elems_old == 1) {
    num_children = t8dg_global_values_get_num_children (t8dg_values_get_global_values (adapt_data->dg_values, itree, first_ielem_old));
    T8DG_ASSERT (num_children == num_elems_new);
    for (ichild = 0; ichild < num_children; ichild++) {
      t8dg_values_set_element_adapt (adapt_data->dg_values, itree, first_ielem_new + ichild);
      t8dg_values_transform_parent_dof_to_child_dof (adapt_data->dg_values, adapt_data->dof_values, adapt_data->dof_values_adapt, itree,
                                                     first_ielem_old, first_ielem_new + ichild, ichild);
    }
  }
  else {
    num_children =
      t8dg_global_values_get_num_children (t8dg_values_get_global_values_adapt (adapt_data->dg_values, itree, first_ielem_new));
    T8DG_ASSERT (num_children == num_elems_old && num_elems_new == 1);
    /* Needed before! we calculate the dof_values */
    t8dg_values_set_element_adapt (adapt_data->dg_values, itree, first_ielem_new);
    t8dg_values_transform_child_dof_to_parent_dof (adapt_data->dg_values, adapt_data->dof_values, adapt_data->dof_values_adapt, itree,
                                                   num_children, first_ielem_old, first_ielem_new);
  }
}

/* Coarsens the mesh in order to apply a multigrid preconditioner */
/* This functions coarsens every family of elements on level back down */
int
t8dg_adapt_multigrid_coarsen_finest_level (t8_forest_t forest,
                                           t8_forest_t forest_from,
                                           t8_locidx_t itree, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements,
                                           t8_element_t * elements[])
{
  t8dg_adapt_data_t  *adapt_data;
  adapt_data = (t8dg_adapt_data_t *) t8_forest_get_user_data (forest);

  if (num_elements == 1) {
    /* If it is not a family, nothing happens */
    return 0;
  }
  /* Right now: Refine at max. to the coarsest geometry possible */
  if (ts->t8_element_level (elements[0]) > 0) {
    return -1;
  }

  /* In case the element should not be coarsened and stays the same */
  return 0;
}
