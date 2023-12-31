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
#include "t8dg_geometry.h"
#include <t8_forest.h>
#include <t8_eclass.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_vec.h>

static void
t8dg_geometry_transform_reference_vertex_to_coarse_vertex (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest,
                                                           const t8_gloidx_t iglobaltree, const t8_element_t * element,
                                                           const double reference_vertex[DIM3], double coarse_vertex[DIM3])
{
  int                 level;
  double              translation_vector[3] = { 0, 0, 0 };
  double              scaling_factor;

  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;

  eclass = t8dg_eclass_from_gloidx_element (forest, iglobaltree, element);
  scheme = t8_forest_get_eclass_scheme (forest, eclass);

  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -level);

  scheme->t8_element_vertex_reference_coords (element, 0, translation_vector);

  /* For triangle reflection about x=y also needed */

  t8_vec_axpyz (reference_vertex, translation_vector, coarse_vertex, scaling_factor);
  /*hx+x_0 */
}

void
t8dg_geometry_transform_reference_vertex_to_image_vertex (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest,
                                                          const t8_gloidx_t iglobaltree, const t8_element_t * element,
                                                          const double reference_vertex[3], double image_vertex[3])
{
  double              coarse_vertex[3];
  /* transform into coarse reference element */
  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (coarse_geometry, forest, iglobaltree, element, reference_vertex,
                                                             coarse_vertex);

  /* tree vertices are application geometry for linear geometry */
  t8dg_coarse_geometry_apply (coarse_geometry, forest, iglobaltree, coarse_vertex, image_vertex);
}

double
t8dg_geometry_calculate_sqrt_gram_determinant (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest,
                                               const t8_gloidx_t iglobaltree, const t8_element_t * element,
                                               const double reference_vertex[3])
{
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;

  int                 dim;
  int                 level;
  double              scaling_factor;
  double              coarse_vertex[3];

  eclass = t8dg_eclass_from_gloidx_element (forest, iglobaltree, element);
  scheme = t8_forest_get_eclass_scheme (forest, eclass);

  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (coarse_geometry, forest, iglobaltree, element, reference_vertex,
                                                             coarse_vertex);

  dim = t8_eclass_to_dimension[eclass];
  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -dim * level);

  return scaling_factor * t8dg_coarse_geometry_calculate_sqrt_gram_determinant (coarse_geometry, forest, iglobaltree, coarse_vertex);
}

void
t8dg_geometry_calculate_transformed_gradient_tangential_vector (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest,
                                                                const t8_gloidx_t iglobaltree, const t8_element_t * element,
                                                                const double reference_vertex[3],
                                                                const double reference_tangential_vector[3],
                                                                double transformed_gradient_tangential_vector[3])
{
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;

  int                 level;
  double              scaling_factor;
  double              coarse_vertex[3];
  double              coarse_tangential_vector[3];

  eclass = t8dg_eclass_from_gloidx_element (forest, iglobaltree, element);
  scheme = t8_forest_get_eclass_scheme (forest, eclass);

  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (coarse_geometry, forest, iglobaltree, element, reference_vertex,
                                                             coarse_vertex);

  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -level);

  t8_vec_axb (reference_tangential_vector, coarse_tangential_vector, 1. / scaling_factor, 0);

  t8dg_coarse_geometry_calculate_gradient_tangential_vector (coarse_geometry, forest,
                                                             iglobaltree, coarse_vertex, coarse_tangential_vector,
                                                             transformed_gradient_tangential_vector);
}

void
t8dg_geometry_calculate_normal_vector (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest, const t8_gloidx_t iglobaltree,
                                       const t8_element_t * element, const int iface, const double reference_vertex[3],
                                       double image_normal_vector[3])
{
  t8_eclass_t         eclass;
  double              coarse_vertex[3] = { 0, 0, 0 };
  double              coarse_normal_vector[3] = { 0, 0, 0 };
  int                 side;
  int                 normal_direction;

  eclass = t8dg_eclass_from_gloidx_element (forest, iglobaltree, element);

  if (eclass == T8_ECLASS_LINE || eclass == T8_ECLASS_QUAD || eclass == T8_ECLASS_HEX) {
    /* For Hypercubes, there is no change of the normal vector from reference element to coarse element */
    side = iface % 2;
    normal_direction = iface / 2;
    coarse_normal_vector[normal_direction] = (side) ? 1 : -1;
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (coarse_geometry, forest, iglobaltree, element, reference_vertex,
                                                             coarse_vertex);

  t8dg_coarse_geometry_transform_normal_vector (coarse_geometry, forest, iglobaltree,
                                                coarse_vertex, coarse_normal_vector, image_normal_vector);
  t8_vec_ax (image_normal_vector, 1. / t8_vec_norm (image_normal_vector));
}

double
t8dg_geometry_calculate_face_sqrt_gram_determinant (const t8dg_coarse_geometry_t * coarse_geometry, t8_forest_t forest,
                                                    const t8_gloidx_t iglobaltree, const t8_element_t * element, const int iface,
                                                    const double reference_vertex[3])
{
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *scheme;

  int                 dim;
  int                 level;
  double              scaling_factor;
  double              coarse_vertex[3];

  eclass = t8dg_eclass_from_gloidx_element (forest, iglobaltree, element);
  scheme = t8_forest_get_eclass_scheme (forest, eclass);

  t8dg_geometry_transform_reference_vertex_to_coarse_vertex (coarse_geometry, forest, iglobaltree, element, reference_vertex,
                                                             coarse_vertex);

  dim = t8_eclass_to_dimension[eclass];
  level = scheme->t8_element_level (element);
  scaling_factor = pow (2, -(dim - 1) * level);

  return scaling_factor * t8dg_coarse_geometry_calculate_sqrt_face_gram_determinant (coarse_geometry, forest, iglobaltree, iface,
                                                                                     coarse_vertex);
}
