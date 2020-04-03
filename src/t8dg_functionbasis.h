/*
 * t8dg_functionbasis.h
 *
 *  Created on: Apr 2, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_FUNCTIONBASIS_H_
#define SRC_T8DG_FUNCTIONBASIS_H_

#include "t8dg_vertexset.h"
#include <sc_containers.h>

typedef enum t8dg_functionbasis_type
{
  T8DG_LAGRANGE_LGL,
  T8DG_LAGRANGE_GL
} t8dg_functionbasis_type_t;

typedef struct t8dg_functionbasis t8dg_functionbasis_t;

t8dg_functionbasis_t *t8dg_functionbasis_new_Lagrange (t8dg_vertexset_t * vertexset);

void                t8dg_functionbasis_destroy (t8dg_functionbasis_t ** pfunctionbasis);

void                t8dg_functionbasis_apply_derivative_matrix_transpose (sc_array_t * dof_values,
                                                                          sc_array_t * derivative_dof_values,
                                                                          const t8dg_functionbasis_t * functionbasis);

t8dg_functionbasis_type_t t8dg_functionbasis_get_type (const t8dg_functionbasis_t * functionbasis);

void                t8dg_functionbasis_get_vertex (double vertex[3], const t8dg_functionbasis_t * functionbasis, const int idof);

int                 t8dg_functionbasis_get_num_dof (const t8dg_functionbasis_t * functionbasis);

int                 t8dg_functionbasis_get_dim (const t8dg_functionbasis_t * functionbasis);

#endif /* SRC_T8DG_FUNCTIONBASIS_H_ */
