/*
 * t8dg_local_precomputed_values.h
 *
 *  Created on: Apr 5, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_LOCAL_PRECOMPUTED_VALUES_H_
#define SRC_T8DG_LOCAL_PRECOMPUTED_VALUES_H_

#include <t8.h>
#include <t8_forest.h>
#include <t8_element.h>

#include "t8dg_quadrature.h"

typedef struct t8dg_local_precomputed_values t8dg_local_precomputed_values_t;

void                t8dg_local_precomputed_values_copy_element_values
  (t8dg_local_precomputed_values_t * incoming_values, t8_locidx_t incoming_idata,
   t8dg_local_precomputed_values_t * outgoing_values, t8_locidx_t outgoing_idata);

void                t8dg_local_precomputed_values_set_element (t8dg_local_precomputed_values_t * values,
                                                               t8_forest_t forest, t8_locidx_t itree, t8_eclass_scheme_c * scheme,
                                                               t8_locidx_t ielement, t8dg_quadrature_t * quadrature);

t8dg_local_precomputed_values_t *t8dg_local_precomputed_values_new (t8dg_quadrature_t * quadrature, t8_locidx_t num_local_elems);

void                t8dg_local_precomputed_values_destroy (t8dg_local_precomputed_values_t ** pvalues);

/** Apply the transformation from fine reference element to the subset of the coarse reference element to a vertex*/
void                t8dg_local_precomputed_values_fine_to_coarse_geometry (const double refined_element_vertex[DIM3],
                                                                           double coarse_element_vertex[DIM3],
                                                                           t8_eclass_scheme_c * scheme, const t8_element_t * element);
/*TODO: eclass scheme as const not possible since 'passing as ‘this’ argument discards qualifiers' */

double             *t8dg_local_precomputed_values_get_transformed_gradient_tangential_vector (t8dg_local_precomputed_values_t * values,
                                                                                              t8_locidx_t idata, t8dg_quad_idx_t iquad,
                                                                                              int idim);

double             *t8dg_local_precomputed_values_get_face_normal_vector (t8dg_local_precomputed_values_t * values, t8_locidx_t idata,
                                                                          int iface, t8dg_quad_idx_t iquad);

void                t8dg_local_precomputed_values_element_multiply_trafo_quad_weight (const t8dg_local_precomputed_values_t * local_values,
                                                                                      sc_array_t * array, t8_locidx_t idata);

void                t8dg_local_precomputed_values_element_divide_trafo_quad_weight (const t8dg_local_precomputed_values_t * local_values,
                                                                                    sc_array_t * array, t8_locidx_t idata);

void                t8dg_local_precomputed_values_partition (t8_forest_t forest_old, t8_forest_t forest_partition,
                                                             t8dg_local_precomputed_values_t * local_values_old,
                                                             t8dg_local_precomputed_values_t * local_values_partition);

#endif /* SRC_T8DG_LOCAL_PRECOMPUTED_VALUES_H_ */