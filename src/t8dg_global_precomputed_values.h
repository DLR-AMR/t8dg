/*
 * t8dg_global_precomputed_values.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */

#include "t8dg_quadrature.h"
#include "t8dg_functionbasis.h"

#ifndef SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_
#define SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_

typedef struct t8dg_global_precomputed_values t8dg_global_precomputed_values_t;

t8dg_global_precomputed_values_t *t8dg_global_precomputed_values_new_1D_LGL (const int number_of_LGL_vertices);

void                t8dg_global_precomputed_values_destroy (t8dg_global_precomputed_values_t ** pvalues);

void                t8dg_global_precomputed_values_transform_element_dof_to_face_quad (const t8dg_global_precomputed_values_t * values,
                                                                                       const int iface,
                                                                                       const sc_array_t * element_dof_array,
                                                                                       sc_array_t * face_quad_array);

void                t8dg_global_precomputed_values_transform_face_quad_to_element_dof (t8dg_global_precomputed_values_t * values,
                                                                                       const int iface,
                                                                                       const sc_array_t * face_quad_array,
                                                                                       sc_array_t * element_dof_array);

t8dg_quadrature_t  *t8dg_global_precomputed_values_get_quadrature (t8dg_global_precomputed_values_t * values);
t8dg_functionbasis_t *t8dg_global_precomputed_values_get_functionbasis (t8dg_global_precomputed_values_t * values);

int                 t8dg_global_precomputed_values_get_num_faces (t8dg_global_precomputed_values_t * values);

#endif /* SRC_T8DG_GLOBAL_PRECOMPUTED_VALUES_H_ */