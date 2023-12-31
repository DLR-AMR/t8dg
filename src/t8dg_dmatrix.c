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
#include "t8dg_dmatrix.h"
#include <sc_containers.h>

#ifdef T8DG_MKL_BLAS
#include <mkl.h>
#else
#include <cblas.h>
#endif

struct t8dg_dmatrix
{
  double             *values;
  int                 nrows;
  int                 ncolumns;
  int                 owner;
};

double
t8dg_dmatrix_at (const t8dg_dmatrix_t * matrix, const int irow, const int icolumn)
{
  T8DG_ASSERT (irow >= 0 && irow < matrix->nrows);
  T8DG_ASSERT (icolumn >= 0 && icolumn < matrix->ncolumns);
  return matrix->values[irow * matrix->ncolumns + icolumn];
}

void
t8dg_dmatrix_set_at (t8dg_dmatrix_t * matrix, const int irow, const int icolumn, const double value)
{
  T8DG_ASSERT (irow >= 0 && irow < matrix->nrows);
  T8DG_ASSERT (icolumn >= 0 && icolumn < matrix->ncolumns);
  matrix->values[irow * matrix->ncolumns + icolumn] = value;
}

void
t8dg_dmatrix_mult_sc_array (const t8dg_dmatrix_t * A, const sc_array_t * x, sc_array_t * b)
{
  /*TODO ASSERT */
  cblas_dgemv (CblasRowMajor, CblasNoTrans, A->nrows, A->ncolumns, 1, A->values, A->ncolumns, (double *) x->array, 1, 0,
               (double *) b->array, 1);
}

void
t8dg_dmatrix_transpose_mult_sc_array (const t8dg_dmatrix_t * A, const sc_array_t * x, sc_array_t * b)
{
  /*TODO ASSERT */
  cblas_dgemv (CblasRowMajor, CblasTrans, A->nrows, A->ncolumns, 1, A->values, A->ncolumns, (double *) x->array, 1, 0, (double *) b->array,
               1);
}

void
t8dg_dmatrix_scale_row (t8dg_dmatrix_t * matrix, const int irow, const double alpha)
{
  T8DG_ASSERT (irow >= 0 && irow < matrix->nrows);
  int                 icolumn;
  for (icolumn = 0; icolumn < matrix->ncolumns; icolumn++) {
    matrix->values[irow * matrix->ncolumns + icolumn] *= alpha;
  }
}

t8dg_dmatrix_t     *
t8dg_dmatrix_new (int nrows, int ncolumns)
{
  T8DG_ASSERT (nrows > 0 && ncolumns > 0);
  t8dg_dmatrix_t     *matrix = T8DG_ALLOC (t8dg_dmatrix_t, 1);
  matrix->values = T8DG_ALLOC (double, nrows * ncolumns);
  matrix->nrows = nrows;
  matrix->ncolumns = ncolumns;
  matrix->owner = 1;
  return matrix;
}

t8dg_dmatrix_t     *
t8dg_dmatrix_new_zero (int nrows, int ncolumns)
{
  T8DG_ASSERT (nrows > 0 && ncolumns > 0);
  t8dg_dmatrix_t     *matrix = T8DG_ALLOC (t8dg_dmatrix_t, 1);
  matrix->values = T8DG_ALLOC_ZERO (double, nrows * ncolumns);
  matrix->nrows = nrows;
  matrix->ncolumns = ncolumns;
  matrix->owner = 1;
  return matrix;
}

t8dg_dmatrix_t     *
t8dg_dmatrix_new_data (int nrows, int ncolumns, double *data)
{
  T8DG_ASSERT (nrows > 0 && ncolumns > 0);
  t8dg_dmatrix_t     *matrix = T8DG_ALLOC (t8dg_dmatrix_t, 1);
  matrix->values = data;
  matrix->nrows = nrows;
  matrix->ncolumns = ncolumns;
  matrix->owner = 0;
  return matrix;
}

void
t8dg_dmatrix_destroy (t8dg_dmatrix_t ** pmatrix)
{
  T8DG_ASSERT (pmatrix != NULL);
  t8dg_dmatrix_t     *matrix = *pmatrix;
  T8DG_ASSERT (matrix != NULL);

  if (matrix->owner) {
    T8DG_FREE (matrix->values);
  }
  matrix->values = NULL;

  matrix->ncolumns = -1;
  matrix->nrows = -1;

  T8DG_FREE (matrix);
  *pmatrix = NULL;
}

void
t8dg_square3D_matrix_invert_sub_matrix (t8dg_square_3D_matrix_t matrix_invers, t8dg_square_3D_matrix_t matrix, int dim)
{
  T8_ASSERT (dim > 0 && dim <= DIM3);
  double              det;
  t8dg_square3D_matrix_determinant_sub_matrix (&det, matrix, dim);
  if (dim == 1) {
    matrix_invers[0][0] = 1. / det;
  }
  else if (dim == 2) {
    matrix_invers[0][0] = 1. / det * matrix[1][1];
    matrix_invers[1][1] = 1. / det * matrix[0][0];
    matrix_invers[0][1] = -1. / det * matrix[0][1];
    matrix_invers[1][0] = -1. / det * matrix[1][0];
  }
  else if (dim == 3) {
    SC_ABORT ("not yet implemented");
  }
}

void
t8dg_square3D_matrix_determinant_sub_matrix (double *det, t8dg_square_3D_matrix_t matrix, int dim)
{
  T8_ASSERT (dim > 0 && dim <= DIM3);
  if (dim == 1) {
    *det = matrix[0][0];
  }
  else if (dim == 2) {
    *det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  }
  else if (dim == 3) {
    SC_ABORT ("not yet implemented");
  }
}

void
t8dg_square3D_matrix_copy (t8dg_square_3D_matrix_t matrix_result, t8dg_square_3D_matrix_t matrix, int dim)
{
  int                 ixdim, iydim;
  for (ixdim = 0; ixdim < dim; ixdim++) {
    for (iydim = 0; iydim < dim; iydim++) {
      matrix_result[ixdim][iydim] = matrix[ixdim][iydim];
    }
  }
}

void
t8dg_square3D_matrix_scale (t8dg_square_3D_matrix_t matrix, double scaling_factor, int dim)
{
  int                 ixdim, iydim;
  for (ixdim = 0; ixdim < dim; ixdim++) {
    for (iydim = 0; iydim < dim; iydim++) {
      matrix[ixdim][iydim] *= scaling_factor;
    }
  }
}

void
t8dg_square3D_matrix_apply_rotation_reflection_matrix_to_matrix (t8dg_square_3D_matrix_t matrix_result,
                                                                 t8dg_square_3D_matrix_t matrix, int idx_rotation_reflection, int dim)
{
  if (idx_rotation_reflection == 0) {
    t8dg_square3D_matrix_copy (matrix_result, matrix, dim);
  }
}
