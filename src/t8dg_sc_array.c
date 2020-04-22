/*
 * t8dg_sc_array.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include <t8.h>
#include "t8dg.h"
#include <sc_containers.h>

void
t8dg_sc_array_copy_only_at_indices (sc_array_t * incoming_array, t8_locidx_t incoming_idata,
                                    sc_array_t * outgoing_array, t8_locidx_t outgoing_idata)
{
  T8DG_ASSERT (incoming_array->elem_size == outgoing_array->elem_size);
  memcpy (t8_sc_array_index_locidx (outgoing_array, outgoing_idata),
          t8_sc_array_index_locidx (incoming_array, incoming_idata), incoming_array->elem_size);
}

sc_array_t         *
t8dg_sc_array_duplicate (const sc_array_t * src)
{
  sc_array_t         *dest;
  dest = sc_array_new_count (src->elem_size, src->elem_count);
  return dest;
}

sc_array_t         *
t8dg_sc_array_block_double_new_view (sc_array_t * src, t8_locidx_t idata)
{
  T8DG_ASSERT (src->elem_size > 0 && src->elem_size % 8 == 0);
  return sc_array_new_data (t8_sc_array_index_locidx (src, idata), sizeof (double), src->elem_size / sizeof (double));
}

sc_array_t         *
t8dg_sc_array_clone (sc_array_t * src)
{
  sc_array_t         *dest = sc_array_new_count (src->elem_size, src->elem_count);
  sc_array_copy (dest, src);
  return dest;
}

void
t8dg_sc_array_copy (const sc_array_t * src, sc_array_t * dest)
{
  SC_ASSERT (dest->elem_size == src->elem_size);
  SC_ASSERT (dest->elem_count == src->elem_count);

  memcpy (dest->array, src->array, src->elem_count * src->elem_size);
}

void
t8dg_sc_array_block_double_axpy (const double a, const sc_array_t * x, sc_array_t * y)
{
  T8DG_ASSERT (x != NULL && y != NULL);
  T8DG_ASSERT (x->array != NULL && y->array != NULL);
  T8DG_ASSERT (x->elem_size % sizeof (double) == 0 && y->elem_size % sizeof (double) == 0);
  T8DG_ASSERT (x->elem_count * x->elem_size == y->elem_count * y->elem_size);
  double             *x_double, *y_double;
  int                 double_count, i;

  /*View array as double array */
  x_double = (double *) x->array;
  y_double = (double *) y->array;
  double_count = x->elem_size / sizeof (double) * x->elem_count;        /*total number of doubles */

  for (i = 0; i < double_count; i++) {
    y_double[i] = a * x_double[i] + y_double[i];
  }
}

void
t8dg_sc_array_block_double_axpyz (double a, const sc_array_t * x, const sc_array_t * y, sc_array_t * z)
{
  T8DG_ASSERT (z != NULL && x != NULL && y != NULL);
  T8DG_ASSERT (z->array != NULL && x->array != NULL && y->array != NULL);
  T8DG_ASSERT (x->elem_size % sizeof (double) == 0 && y->elem_size % sizeof (double) == 0 && z->elem_size % sizeof (double) == 0);
  T8DG_ASSERT (x->elem_count == y->elem_count && x->elem_count == z->elem_count);
  T8DG_ASSERT (x->elem_size == y->elem_size && x->elem_size == z->elem_size);

  /*View array as double array */
  double             *x_double, *y_double, *z_double;
  size_t              double_count, idouble;

  x_double = (double *) (x->array);
  y_double = (double *) (y->array);
  z_double = (double *) (z->array);

  double_count = (x->elem_size / sizeof (double)) * x->elem_count;      /*total number of doubles */

  for (idouble = 0; idouble < double_count; idouble++) {
    z_double[idouble] = a * x_double[idouble] + y_double[idouble];
  }
}

void
t8dg_sc_array_swap (sc_array_t ** parray1, sc_array_t ** parray2)
{
  sc_array_t         *temp;
  temp = *parray1;
  *parray1 = *parray2;
  *parray2 = temp;
}

void
t8dg_sc_array_block_double_debug_print (sc_array_t * array)
{
#ifdef T8_ENABLE_DEBUG
  size_t              irow, icolumn;
  double             *row_array;
  for (irow = 0; irow < array->elem_count; irow++) {
    row_array = (double *) t8_sc_array_index_locidx (array, irow);
    for (icolumn = 0; icolumn < array->elem_size / 8; icolumn++) {
      printf ("%f  ,  ", row_array[icolumn]);
    }
    printf ("\n");
  }
  printf ("\n");
#endif
}

void
t8dg_sc_array_block_double_print (sc_array_t * array)
{
  size_t              irow, icolumn;
  double             *row_array;
  for (irow = 0; irow < array->elem_count; irow++) {
    row_array = (double *) t8_sc_array_index_locidx (array, irow);
    for (icolumn = 0; icolumn < array->elem_size / 8; icolumn++) {
      printf ("%f  ,  ", row_array[icolumn]);
    }
    printf ("\n");
  }
  printf ("\n");
}

void
t8dg_sc_array_block_double_set_zero (sc_array_t * array)
{
  T8DG_ASSERT (array->elem_size % sizeof (double) == 0);
  size_t              number_of_doubles = array->elem_size / sizeof (double) * array->elem_count;
  size_t              idouble;
  double             *array_double;
  array_double = (double *) array->array;
  for (idouble = 0; idouble < number_of_doubles; idouble++) {
    array_double[idouble] = 0;
  }
}

void
t8dg_sc_array_block_square_values (sc_array_t * src, sc_array_t * dest)
{
  T8DG_ASSERT (src != NULL && dest != NULL);
  T8DG_ASSERT (src->array != NULL && dest->array != NULL);
  T8DG_ASSERT (src->elem_size % sizeof (double) == 0 && dest->elem_size % sizeof (double) == 0);
  T8DG_ASSERT (src->elem_count == dest->elem_count);
  T8DG_ASSERT (src->elem_size == dest->elem_size);
  double             *src_array_double;
  double             *dest_array_double;
  size_t              idx;
  size_t              double_count;
  double_count = src->elem_count * src->elem_size / sizeof (double);
  src_array_double = (double *) sc_array_index (src, 0);
  dest_array_double = (double *) sc_array_index (dest, 0);
  for (idx = 0; idx < double_count; idx++) {
    dest_array_double[idx] = SC_SQR (src_array_double[idx]);
  }
}
