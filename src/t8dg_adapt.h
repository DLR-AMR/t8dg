/*
 *  This file is part of t8dg.
 *  t8dg is a discontinuous galerkin solver on adaptive meshes
 *  that uses the adaptive meshing library t8code to distribute
 *  the workload equally among parallel processes to achieve good scaling.
 *
 *  Copyright (C) 2020 the developers
 *
 *  t8dg is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  t8dg is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with t8dg.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SRC_T8DG_ADAPT_H_
#define SRC_T8DG_ADAPT_H_

#include "t8dg.h"
#include "t8dg_dof.h"
#include "t8dg_values.h"

T8DG_EXTERN_C_BEGIN ();

typedef struct t8dg_adapt_data
{
  int                 initial_refinement_level;
  int                 maximum_refinement_level;
  int                 minimum_refinement_level;

  int                 adapt_freq;

  t8dg_dof_values_t  *dof_values;
  t8dg_dof_values_t  *dof_values_adapt;

  t8dg_values_t      *dg_values;

  t8_forest_adapt_t   adapt_fn;

  t8dg_scalar_function_3d_time_fn source_sink_fn;
  void               *source_sink_data;
  t8dg_dof_values_t  *source_sink_dof;
  double              time;
  int                 dim;
} t8dg_adapt_data_t;

t8_forest_adapt_t   t8dg_adapt_fn_arg (int adapt_arg);

t8dg_adapt_data_t  *t8dg_adapt_data_new (t8dg_values_t * dg_values, int initial_level, int min_level, int max_level, int adapt_fn_arg,
                                         int adapt_freq, t8dg_scalar_function_3d_time_fn source_sink_fn, void *source_sink_data);

void                t8dg_adapt_data_set_time (t8dg_adapt_data_t * adapt_data, double time);

void                t8dg_adapt_data_interpolate_source_fn (t8dg_adapt_data_t * adapt_data);

void                t8dg_adapt_data_destroy (t8dg_adapt_data_t ** p_adapt_data);

void                t8dg_adapt_replace
  (t8_forest_t forest_old, t8_forest_t forest_new,
   t8_locidx_t itree,
   t8_eclass_scheme_c * ts, int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new);

int                 t8dg_adapt_rel_min_max
  (t8_forest_t forest, t8_forest_t forest_from,
   t8_locidx_t itree, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

int
 
 
 
 
 
 
 
 
t8dg_adapt_mass (t8_forest_t forest,
                 t8_forest_t forest_from,
                 t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

int
 
 
 
 
 
 
 
 
t8dg_adapt_smooth_indicator (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

int
 
 
 
 
 
 
 
 
t8dg_adapt_smooth_indicator_hypercube (t8_forest_t forest,
                                       t8_forest_t forest_from,
                                       t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements,
                                       t8_element_t * elements[]);

int
 
 
 
 
 
 
 
 
t8dg_adapt_uniform (t8_forest_t forest,
                    t8_forest_t forest_from,
                    t8_locidx_t itree, t8_locidx_t ielement, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

T8DG_EXTERN_C_END ();

#endif
