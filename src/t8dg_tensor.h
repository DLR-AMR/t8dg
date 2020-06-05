/*
 * t8dg_tensor.h
 *
 *  Created on: Apr 24, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_TENSOR_H_
#define SRC_T8DG_TENSOR_H_

#include <sc_containers.h>
#include "t8dg.h"

T8DG_EXTERN_C_BEGIN ();

void                t8dg_tensor_array_extract_vector (sc_array_t * tensor_array, const int ivector, const int stride, sc_array_t * vector);

void                t8dg_tensor_array_inject_vector (sc_array_t * vector, const int ivector, const int stride, sc_array_t * tensor_array);

/*deprecated?*/
int                 t8dg_tensor_mult_other_lengths (const int num_tensor, const int tensor_length[DIM3], const int itensor);

void                t8dg_tensor_transform_tensoridx (const int idx, const int tensordims[DIM3], int itensor[DIM3]);     /*MACRO? */

t8_eclass_t         t8dg_tensor_eclass (t8_eclass_t eclass_tensor1, t8_eclass_t eclass_tensor2);

T8DG_EXTERN_C_END ();

#endif /* SRC_T8DG_TENSOR_H_ */
