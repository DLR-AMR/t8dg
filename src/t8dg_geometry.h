/*
 * geometry.h
 *
 *  Created on: Apr 27, 2020
 *      Author: lukas
 */

#ifndef SRC_T8DG_GEOMETRY_H_
#define SRC_T8DG_GEOMETRY_H_

#include <t8.h>
#include <t8_forest.h>
#include "t8dg_coarse_geometry.h"

typedef struct t8dg_geometry_transformation_data
{
  const t8dg_coarse_geometry_t *coarse_geometry;
  const t8_forest_t   forest;
  t8_locidx_t         itree;
  t8_locidx_t         ielement;

} t8dg_geometry_transformation_data_t;

void                t8dg_geometry_transform_reference_vertex_to_image_vertex
  (const t8dg_geometry_transformation_data_t * geometry_data, const double reference_vertex[3], double image_vertex[3]);

double              t8dg_geometry_calculate_gram_determinant
  (const t8dg_geometry_transformation_data_t * geometry_data, const double reference_vertex[3]);

void                t8dg_geometry_calculate_transformed_gradient_tangential_vector
  (const t8dg_geometry_transformation_data_t * geometry_data,
   const double reference_vertex[3], const double reference_tangential_vector[3], double transformed_gradient_tangential_vector[3]);

void                t8dg_geometry_calculate_normal_vector
  (const t8dg_geometry_transformation_data_t * geometry_data, const int iface,
   const double reference_vertex[3], double image_normal_vector[3]);

double              t8dg_geometry_calculate_face_gram_determinant
  (const t8dg_geometry_transformation_data_t * geometry_data, const int iface, const double reference_vertex[3]);

#endif /* SRC_T8DG_GEOMETRY_H_ */