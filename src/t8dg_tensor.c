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



#include "t8dg.h"
#include <sc_containers.h>
#include <t8_eclass.h>
#include "t8dg_tensor.h"

t8_eclass_t
t8dg_tensor_eclass (t8_eclass_t eclass_tensor1, t8_eclass_t eclass_tensor2)
{
  if (eclass_tensor1 == T8_ECLASS_VERTEX) {
    return eclass_tensor2;
  }
  if (eclass_tensor2 == T8_ECLASS_VERTEX) {
    return eclass_tensor1;
  }
  if (eclass_tensor1 == T8_ECLASS_LINE && eclass_tensor2 == T8_ECLASS_LINE) {
    return T8_ECLASS_QUAD;
  }
  if ((eclass_tensor1 == T8_ECLASS_LINE && eclass_tensor2 == T8_ECLASS_QUAD) ||
      (eclass_tensor1 == T8_ECLASS_QUAD && eclass_tensor2 == T8_ECLASS_LINE)) {
    return T8_ECLASS_HEX;
  }
  if ((eclass_tensor1 == T8_ECLASS_LINE && eclass_tensor2 == T8_ECLASS_TRIANGLE) ||
      (eclass_tensor1 == T8_ECLASS_TRIANGLE && eclass_tensor2 == T8_ECLASS_LINE)) {
    return T8_ECLASS_PRISM;
  }
  T8DG_ABORT ("Not possible in 3D");
}

void
t8dg_tensor_array_extract_vector (sc_array_t * tensor_array, const int ivector, const int stride, sc_array_t * vector)
{
  T8DG_ASSERT (stride > 0);
  T8DG_ASSERT (ivector >= 0);

  int                 ivec_element, tensor_array_index, vector_length;
  vector_length = vector->elem_count;
  tensor_array_index = (ivector / stride) * vector_length * stride + (ivector % stride);

  for (ivec_element = 0; ivec_element < vector_length; ivec_element++) {
    *(double *) sc_array_index_int (vector, ivec_element) = *(double *) sc_array_index_int (tensor_array, tensor_array_index);
    tensor_array_index += stride;
  }
}

void
t8dg_tensor_transform_tensoridx (int idx, const int tensordims[DIM3], int tensoridx[DIM3])
{
  /*TODO: Check */
  int                 numtensor = 0, itensor;
  while (numtensor < 3 && tensordims[numtensor] > 0) {
    numtensor++;
  }
  for (itensor = 0; itensor < numtensor; itensor++) {
    tensoridx[itensor] = idx % tensordims[itensor];
    idx /= tensordims[itensor];
  }
}

int
t8dg_tensor_mult_other_lengths (const int num_tensor, const int tensor_length[DIM3], const int itensor)
{
  int                 result = 1;
  int                 iothertensor;
  for (iothertensor = 0; iothertensor < num_tensor; iothertensor++) {
    if (iothertensor != itensor) {
      result *= tensor_length[iothertensor];
    }
  }
  return result;
}

void
t8dg_tensor_array_inject_vector (sc_array_t * vector, const int ivector, const int stride, sc_array_t * tensor_array)
{
  T8DG_ASSERT (ivector >= 0);
  T8DG_ASSERT (stride > 0);

  int                 ivec_element, tensor_array_index, vector_length;
  vector_length = vector->elem_count;
  tensor_array_index = (ivector / stride) * vector_length * stride + (ivector % stride);

  for (ivec_element = 0; ivec_element < vector_length; ivec_element++) {
    *(double *) sc_array_index_int (tensor_array, tensor_array_index) = *(double *) sc_array_index_int (vector, ivec_element);
    tensor_array_index += stride;
  }
}
