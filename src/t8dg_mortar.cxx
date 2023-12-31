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



#include <t8_forest.h>
#include <sc_containers.h>
#include <t8_element_cxx.hxx>
#include <t8_vec.h>

#include "t8dg_mortar.h"
#include "t8dg.h"
#include "t8dg_global_values.h"
#include "t8dg_values.h"
#include "t8dg_flux.h"
#include "t8dg_flux_implementation.h"
#include "t8dg_sc_array.h"
#include "t8dg_geometry.h"

/** struct used to save the calculated numerical fluxes at quadrature points*/
struct t8dg_mortar
{
  t8_gloidx_t         iglobaltree_minus;
  t8_gloidx_t         iglobaltree_plus;
  t8_element_t       *element_minus;
  t8_element_t       *element_plus[MAX_SUBFACES];

  t8_eclass_t         eclass_minus;
  t8_eclass_t         eclass_plus[MAX_SUBFACES];

  t8_locidx_t         idata_minus;                   /**< Local index of the element corresponding to u_minus */
  t8_locidx_t         idata_plus[MAX_SUBFACES];                    /**< Local index of the element corresponding to u_plus */

  int                 iface_minus;
  int                 iface_plus[MAX_SUBFACES];
  int                 num_subfaces;

  t8dg_functionbasis_t *face_functionbasis;    /**< face functionbasis with orientation from elem_minus*/
  t8dg_functionbasis_t *functionbasis_minus;

  int                 number_face_dof;

  int                 bigface_is_ghost;
  int                 subface_is_local[MAX_SUBFACES];

  t8_eclass_t         eclass_face;
  int                 orientation;
  /*one value for each quadrature point, orientation of u_ */
  t8dg_face_quad_values_t *fluxvalue_minus;                                   /**< value of (cu)*.n at face dof */
  /*orientation and sign of right side */
  t8dg_face_quad_values_t *fluxvalue_plus[MAX_SUBFACES];
  int                 valid;                            /**< indicates wether the fluxes are already newly calculated this timestep*/
  int                 is_boundary;

  t8dg_mortar_array_t *mortar_array;
};

struct t8dg_mortar_array
{
  t8_forest_t         forest;
  sc_array_t         *mortars[MAX_FACES];

  int                 num_local_elements;
  int                 num_total_elements;       /*including ghosts */
  int                 max_num_faces;
  /*for hybrid the number of faces for each tree has to be known */
  t8dg_local_values_t *local_values;
  double              ghost_exchange_time;
};

void
t8dg_mortar_destroy (t8dg_mortar_t ** pmortar)
{
  t8dg_mortar        *mortar = *pmortar;
  int                 isubface;

  if (!mortar->bigface_is_ghost) {
    t8dg_debugf ("destroy fluxvalue_minus %p, idata: %i, iface: %i\n", (void *) mortar->fluxvalue_minus, mortar->idata_minus,
                 mortar->iface_minus);
    sc_array_destroy (mortar->fluxvalue_minus);
    t8dg_debugf ("destroyed fluxvalue_minus, idata: %i, iface: %i\n", mortar->idata_minus, mortar->iface_minus);
  }
  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      sc_array_destroy (mortar->fluxvalue_plus[isubface]);
    }
  }
  t8dg_functionbasis_unref (&mortar->functionbasis_minus);
  t8dg_functionbasis_unref (&mortar->face_functionbasis);

  T8DG_FREE (*pmortar);
  *pmortar = NULL;
}

/*make accessible to ghost elements*/
static int
t8dg_mortar_calculate_face_orientation (const t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface)
{
  t8_element_t       *element;
  t8_eclass_scheme_c *scheme;
  int                 orientation;
  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  t8_locidx_t         icmesh_tree;
  icmesh_tree = t8_forest_ltreeid_to_cmesh_ltreeid (forest, itree);

  scheme = t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest, itree));
  if (scheme->t8_element_is_root_boundary (element, iface)) {
    t8_cmesh_get_face_neighbor (t8_forest_get_cmesh (forest), icmesh_tree, scheme->t8_element_tree_face (element, iface), NULL,
                                &orientation);
    return orientation;
  }
  else {
    return 0;
  }
}

static t8dg_mortar_t *
t8dg_mortar_array_get_mortar (const t8dg_mortar_array_t * mortar_array, const t8_locidx_t idata, const int iface)
{
  return *(t8dg_mortar_t **) sc_array_index (mortar_array->mortars[iface], idata);
}

/*Is only called by local elements*/
static t8dg_mortar_t *
t8dg_mortar_new (t8_forest_t forest, t8_locidx_t itree, t8_locidx_t ielement, int iface, t8dg_mortar_array_t * mortar_array)
{
  t8dg_mortar_t      *mortar;
  t8_locidx_t         idata;
  int                 neigh_itree, neigh_ielement, isubface;
  int                *neigh_ifaces;
  int                 num_neighs;
  int                 own_level;
  t8_element_t       *element, **neigh_elems;
  t8_eclass_scheme_c *neigh_scheme;
  t8_locidx_t        *neigh_idatas;
  t8dg_global_values_t *global_values;
  t8dg_functionbasis_t *element_functionbasis;

  global_values = t8dg_local_values_get_global_values (mortar_array->local_values, itree, ielement);
  element_functionbasis = t8dg_global_values_get_functionbasis (global_values);

  element = t8_forest_get_element_in_tree (forest, itree, ielement);

  own_level = t8_forest_get_eclass_scheme (forest, t8_forest_get_eclass (forest, itree))->t8_element_level (element);

  /* use this function to also get the data-indices of the neighbouring element */
  t8_forest_leaf_face_neighbors (forest, itree, element, &neigh_elems, iface, &neigh_ifaces, &num_neighs, &neigh_idatas, &neigh_scheme, 1);

  if (num_neighs == 1 && (neigh_scheme->t8_element_level (neigh_elems[0]) < own_level)) {
    /*the neighbour element is the bigger one */
    if (neigh_idatas[0] < t8_forest_get_num_element (forest)) {
      /*The neighbour element is local */
      t8_forest_get_element (forest, neigh_idatas[0], &neigh_itree);    /*Only needed for neighbour itree */
      neigh_ielement = neigh_idatas[0] - t8_forest_get_tree_element_offset (forest, neigh_itree);
      mortar = t8dg_mortar_new (forest, neigh_itree, neigh_ielement, neigh_ifaces[0], mortar_array);
    }
    else {
      mortar = t8dg_mortar_array_get_mortar (mortar_array, neigh_idatas[0], neigh_ifaces[0]);
      if (mortar == NULL) {
        /*Create new mortar, since neighbour element is not local */
        mortar = T8DG_ALLOC_ZERO (t8dg_mortar_t, 1);
        mortar->mortar_array = mortar_array;
        mortar->bigface_is_ghost = 1;

        mortar->idata_minus = neigh_idatas[0];
        mortar->element_minus = neigh_elems[0];
        t8_element_t       *element_temp;
        neigh_scheme->t8_element_new (1, &element_temp);
        int                 neigh_iface_temp;
        mortar->iglobaltree_minus =
          t8_forest_element_face_neighbor (forest, itree, element, element_temp, neigh_scheme, iface, &neigh_iface_temp);
        neigh_scheme->t8_element_destroy (1, &element_temp);
        mortar->eclass_minus = t8dg_eclass_from_gloidx_element (forest, mortar->iglobaltree_minus, mortar->element_minus);

        mortar->iface_minus = neigh_ifaces[0];
        mortar->face_functionbasis = t8dg_functionbasis_get_face_functionbasis (element_functionbasis, iface);
        mortar->functionbasis_minus = element_functionbasis;
        mortar->num_subfaces = t8dg_functionbasis_get_num_children (mortar->face_functionbasis);
        t8dg_functionbasis_ref (mortar->face_functionbasis);
        t8dg_functionbasis_ref (mortar->functionbasis_minus);
        mortar->number_face_dof = t8dg_functionbasis_get_num_dof (mortar->face_functionbasis);
        mortar->eclass_face = t8dg_functionbasis_get_eclass (mortar->face_functionbasis);
        mortar->orientation = t8dg_mortar_calculate_face_orientation (forest, itree, ielement, iface);  /*orientation from small to big element */
      }
//      t8dg_debugf ("Adress: %p, idata: %i, num_local: %i\n", (void *) mortar, neigh_idatas[0], t8_forest_get_num_element (forest));
//      t8dg_debugf ("mortar: num_subfaces: %i\n", mortar->num_subfaces);
//      t8dg_debugf ("Adress: %p\n", (void *) mortar);

      t8_eclass_scheme_c *boundary_scheme, *element_scheme;
      t8_element_t       *boundary_element;

      boundary_scheme = t8_forest_get_eclass_scheme (forest, mortar->eclass_face);
      boundary_scheme->t8_element_new (1, &boundary_element);
      element = t8_forest_get_element_in_tree (forest, itree, ielement);
      element_scheme = t8_forest_get_eclass_scheme (forest, t8dg_functionbasis_get_eclass (element_functionbasis));
      element_scheme->t8_element_boundary_face (element, iface, boundary_element, boundary_scheme);

      isubface = boundary_scheme->t8_element_child_id (boundary_element);
      boundary_scheme->t8_element_destroy (1, &boundary_element);

      idata = t8dg_itree_ielement_to_idata (forest, itree, ielement);

//      t8dg_debugf ("Halfmortar set: isubface: %i, idata: %i, iface: %i\n", isubface, idata, iface);

      mortar->idata_plus[isubface] = idata;
      mortar->element_plus[isubface] = t8_forest_get_element_in_tree (forest, itree, ielement);
      mortar->iglobaltree_plus = t8_forest_global_tree_id (forest, itree);
      mortar->eclass_plus[isubface] = t8dg_eclass_from_gloidx_element (forest, mortar->iglobaltree_plus, mortar->element_plus[isubface]);

      mortar->iface_plus[isubface] = iface;
      T8DG_ASSERT (mortar->subface_is_local[isubface] == 0);
      mortar->subface_is_local[isubface] = 1;
      mortar->fluxvalue_plus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
//      t8dg_debugf ("add fluxvalue_plus[isubface] = %p", (void *) mortar->fluxvalue_plus[isubface]);
      mortar->valid = 0;
    }
  }
  else {
    /*own element is bigger or same size */
    mortar = T8DG_ALLOC_ZERO (t8dg_mortar_t, 1);
    if (num_neighs == 0) {
      mortar->is_boundary = 1;
    }
    mortar->mortar_array = mortar_array;
    mortar->bigface_is_ghost = 0;
    idata = t8dg_itree_ielement_to_idata (forest, itree, ielement);

    mortar->iface_minus = iface;
    mortar->idata_minus = idata;
    mortar->element_minus = t8_forest_get_element_in_tree (forest, itree, ielement);
    mortar->iglobaltree_minus = t8_forest_global_tree_id (forest, itree);
    mortar->eclass_minus = t8dg_eclass_from_gloidx_element (forest, mortar->iglobaltree_minus, mortar->element_minus);

    mortar->num_subfaces = num_neighs;

    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      mortar->subface_is_local[isubface] = 1;   /*Since the big face is local, all subfaces are local */
      mortar->idata_plus[isubface] = neigh_idatas[isubface];    /*could be greater than number of local elements -> ghost */
      mortar->iface_plus[isubface] = neigh_ifaces[isubface];    /* get neighbouring face index */

      mortar->element_plus[isubface] = neigh_elems[isubface];

      t8_element_t       *element_temp;
      int                 iface_temp;
      neigh_scheme->t8_element_new (1, &element_temp);
      mortar->iglobaltree_plus = t8_forest_element_face_neighbor (forest, itree, element, element_temp, neigh_scheme, iface, &iface_temp);
      neigh_scheme->t8_element_destroy (1, &element_temp);

      mortar->eclass_plus[isubface] = t8dg_eclass_from_gloidx_element (forest, mortar->iglobaltree_plus, mortar->element_plus[isubface]);
    }
    mortar->face_functionbasis = t8dg_functionbasis_get_face_functionbasis (element_functionbasis, iface);
    mortar->functionbasis_minus = element_functionbasis;
    t8dg_functionbasis_ref (mortar->face_functionbasis);
    t8dg_functionbasis_ref (mortar->functionbasis_minus);
    mortar->number_face_dof = t8dg_functionbasis_get_num_dof (mortar->face_functionbasis);

    mortar->eclass_face = t8dg_functionbasis_get_eclass (mortar->face_functionbasis);
    mortar->orientation = t8dg_mortar_calculate_face_orientation (forest, itree, ielement, iface);      /*orientation from big to small element, TODO: differs */

    /* allocate memory for sc_arrays */
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      mortar->fluxvalue_plus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
    }
    mortar->fluxvalue_minus = sc_array_new_count (sizeof (double), mortar->number_face_dof);
    mortar->valid = 0;
  }
  if (num_neighs > 0) {
    T8_FREE (neigh_elems);
    T8_FREE (neigh_idatas);
    T8_FREE (neigh_ifaces);
  }
  return mortar;
}

static void
t8dg_mortar_calculate_linear_flux3D (t8dg_mortar_t * mortar, t8dg_dof_values_t * dof_values, t8dg_linear_flux3D_fn linear_flux,
                                     void *flux_data, t8dg_numerical_linear_flux3D_fn numerical_flux, void *numerical_flux_data,
                                     double time)
{
  sc_array_t         *elem_dof_values_minus;
  sc_array_t         *elem_dof_values_plus[MAX_SUBFACES];

  sc_array_t         *face_dof_values_minus_big;
  sc_array_t         *face_dof_values_minus[MAX_SUBFACES];
  sc_array_t         *face_dof_values_plus[MAX_SUBFACES];
  sc_array_t         *face_flux_values_minus[MAX_SUBFACES];

  int                 idof;
  int                 isubface;
  double              fluxvalue;
  double              u_minus_val;
  double              u_plus_val;
  double              flux_vec[3];
  double              reference_vertex[3] = { 0, 0, 0 };
  double              image_vertex[3];

  t8dg_mortar_array_t *mortar_array;
  t8dg_functionbasis_t *element_functionbasis;
  element_functionbasis = mortar->functionbasis_minus;
  mortar_array = mortar->mortar_array;

  if (mortar->is_boundary) {
    t8dg_face_dof_values_set_zero (mortar->fluxvalue_minus);
    return;
  }

  face_dof_values_minus_big = sc_array_new_count (sizeof (double), mortar->number_face_dof);
  elem_dof_values_minus = t8dg_dof_values_new_element_dof_values_view_idata_eclass (dof_values, mortar->idata_minus, mortar->eclass_minus);

  t8dg_functionbasis_transform_element_dof_to_face_dof (element_functionbasis, mortar->iface_minus, elem_dof_values_minus,
                                                        face_dof_values_minus_big);

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      elem_dof_values_plus[isubface] =
        t8dg_dof_values_new_element_dof_values_view_idata_eclass (dof_values, mortar->idata_plus[isubface], mortar->eclass_plus[isubface]);

      face_dof_values_plus[isubface] = t8dg_face_dof_values_new (mortar->number_face_dof);
      t8dg_functionbasis_transform_element_dof_to_face_dof (element_functionbasis, mortar->iface_plus[isubface],
                                                            elem_dof_values_plus[isubface], face_dof_values_plus[isubface]);
    }
  }

  if (mortar->num_subfaces > 1) {
    /*Interpolate to children that are local */
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (mortar->subface_is_local[isubface]) {
        face_dof_values_minus[isubface] = t8dg_face_dof_values_new (mortar->number_face_dof);
        t8dg_functionbasis_apply_child_interpolation_matrix (mortar->face_functionbasis, isubface, face_dof_values_minus_big,
                                                             face_dof_values_minus[isubface]);
        t8dg_face_dof_values_orient (face_dof_values_minus[isubface], mortar->eclass_face, mortar->orientation);
      }
    }
  }
  else {
    T8DG_ASSERT (mortar->num_subfaces == 1);
    face_dof_values_minus[0] = face_dof_values_minus_big;
    t8dg_face_dof_values_orient (face_dof_values_minus[0], mortar->eclass_face, mortar->orientation);
  }

  /*Calculate Flux values */
  t8dg_debugf ("calculate flux for mortar %p\n", (void *) mortar);
  if (t8dg_functionbasis_is_lagrange (mortar->face_functionbasis)) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (!mortar->subface_is_local[isubface]) {
        continue;
      }
      t8dg_debugf ("calculate flux for idata %i, iface %i ,isubface %i\n", mortar->idata_plus[isubface], mortar->iface_plus[isubface],
                   isubface);
      face_flux_values_minus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      T8DG_ASSERT (mortar->fluxvalue_plus[isubface]->elem_count == (size_t) mortar->number_face_dof);
      for (idof = 0; idof < mortar->number_face_dof; idof++) {
        double              outward_normal[3] = { 0, 0, 0 };

        t8dg_local_values_get_face_normal_vector_idata_eclass (mortar_array->local_values, mortar->idata_plus[isubface],
                                                               mortar->eclass_plus[isubface], mortar->iface_plus[isubface], idof,
                                                               outward_normal);

        t8_vec_ax (outward_normal, -1); //need outward normal to elem_minus at the face nodal_basis_vertices of elem_plus,

        t8dg_functionbasis_t *functionbasis_plus;
        //t8dg_global_values_t *global_values_plus;
        /*TODO: Wrong, use smaller functionbasis and geometry */
        functionbasis_plus = mortar->functionbasis_minus;

        t8dg_functionbasis_get_lagrange_face_vertex (functionbasis_plus, mortar->iface_plus[isubface], idof, reference_vertex);

        t8dg_geometry_transform_reference_vertex_to_image_vertex (t8dg_local_values_get_coarse_geometry (mortar_array->local_values),
                                                                  mortar_array->forest, mortar->iglobaltree_plus,
                                                                  mortar->element_plus[isubface], reference_vertex, image_vertex);

        linear_flux (image_vertex, flux_vec, time, flux_data);

        u_minus_val = t8dg_face_dof_values_get_value (face_dof_values_minus[isubface], idof);
        u_plus_val = t8dg_face_dof_values_get_value (face_dof_values_plus[isubface], idof);

        fluxvalue = numerical_flux (u_minus_val, u_plus_val, flux_vec, outward_normal, numerical_flux_data);
        t8dg_debugf ("fluxvalue: %f\n", fluxvalue);
        T8DG_ASSERT (fluxvalue == fluxvalue && fabs (fluxvalue) < 1e200);
        *(double *) sc_array_index_int (face_flux_values_minus[isubface], idof) = +fluxvalue;   /*TODO: check! */
        *(double *) sc_array_index_int (mortar->fluxvalue_plus[isubface], idof) = -fluxvalue;
      }
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (face_flux_values_minus[isubface]));
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_plus[isubface]));
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  /*Project on big face if it is local */
  if (!(mortar->bigface_is_ghost)) {
    if (mortar->num_subfaces > 1) {
      /*What is the order of orientation and projection ? */
      t8dg_local_values_transform_orient_face_child_dof_to_parent_dof_hanging_nodes (mortar_array->local_values, face_flux_values_minus,
                                                                                     mortar->fluxvalue_minus, mortar->num_subfaces,
                                                                                     mortar->idata_plus, mortar->eclass_plus,
                                                                                     mortar->iface_plus, mortar->idata_minus,
                                                                                     mortar->eclass_minus, mortar->iface_minus,
                                                                                     mortar->orientation);
    }
    else {
      t8dg_face_dof_values_orient_back (face_flux_values_minus[0], mortar->eclass_face, mortar->orientation);
      sc_array_copy (mortar->fluxvalue_minus, face_flux_values_minus[0]);

    }
  }

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      sc_array_destroy (elem_dof_values_plus[isubface]);
      sc_array_destroy (face_dof_values_minus[isubface]);
      sc_array_destroy (face_flux_values_minus[isubface]);
      sc_array_destroy (face_dof_values_plus[isubface]);
    }
  }
  if (mortar->num_subfaces > 1) {
    sc_array_destroy (face_dof_values_minus_big);

  }
  sc_array_destroy (elem_dof_values_minus);

  mortar->valid = 1;
}

static void
t8dg_mortar_calculate_flux_dof1D (t8dg_mortar_t * mortar, t8dg_dof_values_t * dof_values,
                                  int icomp, t8dg_numerical_flux1D_fn numerical_flux, void *numerical_flux_data, double time)
{
  sc_array_t         *elem_dof_values_minus;
  sc_array_t         *elem_dof_values_plus[MAX_SUBFACES];

  sc_array_t         *face_dof_values_minus_big;
  sc_array_t         *face_dof_values_minus[MAX_SUBFACES];
  sc_array_t         *face_dof_values_plus[MAX_SUBFACES];
  sc_array_t         *face_flux_values_minus[MAX_SUBFACES];

  int                 idof;
  int                 isubface;
  double              fluxvalue;
  double              u_minus_val;
  double              u_plus_val;

  t8dg_mortar_array_t *mortar_array;
  t8dg_functionbasis_t *element_functionbasis = mortar->functionbasis_minus;

  mortar_array = mortar->mortar_array;

  if (mortar->is_boundary) {
    t8dg_face_dof_values_set_zero (mortar->fluxvalue_minus);
    return;
  }

  face_dof_values_minus_big = sc_array_new_count (sizeof (double), mortar->number_face_dof);

  elem_dof_values_minus = t8dg_dof_values_new_element_dof_values_view_idata_eclass (dof_values, mortar->idata_minus, mortar->eclass_minus);

  t8dg_functionbasis_transform_element_dof_to_face_dof (element_functionbasis, mortar->iface_minus, elem_dof_values_minus,
                                                        face_dof_values_minus_big);

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      elem_dof_values_plus[isubface] =
        t8dg_dof_values_new_element_dof_values_view_idata_eclass (dof_values, mortar->idata_plus[isubface], mortar->eclass_plus[isubface]);

      face_dof_values_plus[isubface] = t8dg_face_dof_values_new (mortar->number_face_dof);
      t8dg_functionbasis_transform_element_dof_to_face_dof (element_functionbasis, mortar->iface_plus[isubface],
                                                            elem_dof_values_plus[isubface], face_dof_values_plus[isubface]);
    }
  }

  if (mortar->num_subfaces > 1) {
    /*Interpolate to children that are local */
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (mortar->subface_is_local[isubface]) {
        face_dof_values_minus[isubface] = t8dg_face_dof_values_new (mortar->number_face_dof);
        t8dg_functionbasis_apply_child_interpolation_matrix (mortar->face_functionbasis, isubface, face_dof_values_minus_big,
                                                             face_dof_values_minus[isubface]);
        t8dg_face_dof_values_orient (face_dof_values_minus[isubface], mortar->eclass_face, mortar->orientation);
      }
    }
  }
  else {
    T8DG_ASSERT (mortar->num_subfaces == 1);
    face_dof_values_minus[0] = face_dof_values_minus_big;
    t8dg_face_dof_values_orient (face_dof_values_minus[0], mortar->eclass_face, mortar->orientation);
  }

  /*Calculate Flux values */
  t8dg_debugf ("calculate flux for mortar %p\n", (void *) mortar);
  if (t8dg_functionbasis_is_lagrange (mortar->face_functionbasis)) {
    for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
      if (!mortar->subface_is_local[isubface]) {
        continue;
      }
      t8dg_debugf ("calculate flux for idata %i, iface %i ,isubface %i\n", mortar->idata_plus[isubface], mortar->iface_plus[isubface],
                   isubface);
      face_flux_values_minus[isubface] = sc_array_new_count (sizeof (double), mortar->number_face_dof);
      T8DG_ASSERT (mortar->fluxvalue_plus[isubface]->elem_count == (size_t) mortar->number_face_dof);
      for (idof = 0; idof < mortar->number_face_dof; idof++) {
        double              outward_normal[3] = { 0, 0, 0 };

        t8dg_local_values_get_face_normal_vector_idata_eclass (mortar_array->local_values, mortar->idata_plus[isubface], mortar->eclass_plus[isubface], mortar->iface_plus[isubface], idof, outward_normal);    /*instead of getting all three components and discarding one, write function for single component */

        t8_vec_ax (outward_normal, -1); //need outward normal to elem_minus at the face nodal_basis_vertices of elem_plus,

        u_minus_val = t8dg_face_dof_values_get_value (face_dof_values_minus[isubface], idof);
        u_plus_val = t8dg_face_dof_values_get_value (face_dof_values_plus[isubface], idof);

        int                 reverse_direction;
        /*How to guarantee that left/right is taken correctly? */
        if (numerical_flux == t8dg_numerical_flux1D_left || numerical_flux == t8dg_numerical_flux1D_right) {
          reverse_direction = 1 - (mortar->iface_minus % 2);
          numerical_flux_data = &reverse_direction;
        }
        fluxvalue = numerical_flux (u_minus_val, u_plus_val, outward_normal[icomp], numerical_flux_data);
        t8dg_debugf ("fluxvalue: %f\n", fluxvalue);
        T8DG_ASSERT (fluxvalue == fluxvalue && fabs (fluxvalue) < 1e200);
        *(double *) sc_array_index_int (face_flux_values_minus[isubface], idof) = +fluxvalue;
        *(double *) sc_array_index_int (mortar->fluxvalue_plus[isubface], idof) = -fluxvalue;
      }
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (face_flux_values_minus[isubface]));
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_plus[isubface]));
    }
  }
  else {
    T8DG_ABORT ("Not yet implemented");
  }

  /*Project on big face if it is local */
  if (!(mortar->bigface_is_ghost)) {
    if (mortar->num_subfaces > 1) {
      t8dg_local_values_transform_orient_face_child_dof_to_parent_dof_hanging_nodes (mortar_array->local_values, face_flux_values_minus,
                                                                                     mortar->fluxvalue_minus, mortar->num_subfaces,
                                                                                     mortar->idata_plus, mortar->eclass_plus,
                                                                                     mortar->iface_plus, mortar->idata_minus,
                                                                                     mortar->eclass_minus, mortar->iface_minus,
                                                                                     mortar->orientation);
    }
    else {
      t8dg_face_dof_values_orient_back (face_flux_values_minus[0], mortar->eclass_face, mortar->orientation);
      sc_array_copy (mortar->fluxvalue_minus, face_flux_values_minus[0]);

    }
  }

  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface]) {
      sc_array_destroy (elem_dof_values_plus[isubface]);
      sc_array_destroy (face_dof_values_minus[isubface]);
      sc_array_destroy (face_flux_values_minus[isubface]);
      sc_array_destroy (face_dof_values_plus[isubface]);
    }
  }
  if (mortar->num_subfaces > 1) {
    sc_array_destroy (face_dof_values_minus_big);

  }
  sc_array_destroy (elem_dof_values_minus);

  mortar->valid = 1;
}

static void
t8dg_mortar_array_set_mortar (t8dg_mortar_array_t * mortar_array, const t8_locidx_t idata, const int iface, t8dg_mortar_t * mortar)
{
  T8DG_ASSERT (idata < mortar_array->num_total_elements);
  T8DG_ASSERT (iface < mortar_array->max_num_faces);
  *(t8dg_mortar_t **) sc_array_index (mortar_array->mortars[iface], idata) = mortar;
}

static void
t8dg_mortar_array_set_all_pointers (t8dg_mortar_array_t * mortar_array, t8dg_mortar_t * mortar)
{
  int                 isubface;
  if (!mortar->is_boundary) {
    if (mortar->bigface_is_ghost) {
      for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
        if (mortar->subface_is_local[isubface]) {
          t8dg_debugf ("set mortar %p at: idata:%i, iface:%i\n", (void *) mortar, mortar->idata_plus[isubface],
                       mortar->iface_plus[isubface]);
          t8dg_mortar_array_set_mortar (mortar_array, mortar->idata_plus[isubface], mortar->iface_plus[isubface], mortar);
        }
      }
    }
    else {
      for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
        t8dg_mortar_array_set_mortar (mortar_array, mortar->idata_plus[isubface], mortar->iface_plus[isubface], mortar);
      }
    }
  }
  t8dg_mortar_array_set_mortar (mortar_array, mortar->idata_minus, mortar->iface_minus, mortar);
}

static void
t8dg_mortar_array_set_all_pointers_to_NULL (t8dg_mortar_array_t * mortar_array, t8dg_mortar_t * mortar)
{
  if (!mortar->is_boundary) {
    int                 isubface;
    if (mortar->bigface_is_ghost) {
      for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
        if (mortar->subface_is_local[isubface]) {
          t8dg_mortar_array_set_mortar (mortar_array, mortar->idata_plus[isubface], mortar->iface_plus[isubface], NULL);
        }
      }
    }
    else {
      for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
        t8dg_mortar_array_set_mortar (mortar_array, mortar->idata_plus[isubface], mortar->iface_plus[isubface], NULL);
      }
    }
  }
  t8dg_mortar_array_set_mortar (mortar_array, mortar->idata_minus, mortar->iface_minus, NULL);
}

sc_array_t         *
t8dg_mortar_array_get_oriented_flux (t8dg_mortar_array_t * mortar_array, t8_locidx_t idata, int iface)
{
  t8dg_mortar_t      *mortar;
  int                 isubface;
  t8dg_debugf ("idata: %i, iface: %i \n", idata, iface);
  mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);
  T8DG_ASSERT (mortar != NULL);
  t8dg_debugf ("adress: %p \n", (void *) mortar);
  if (idata == mortar->idata_minus && iface == mortar->iface_minus) {
    T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_minus));
    return mortar->fluxvalue_minus;
  }
  for (isubface = 0; isubface < mortar->num_subfaces; isubface++) {
    if (mortar->subface_is_local[isubface] && mortar->idata_plus[isubface] == idata && mortar->iface_plus[isubface] == iface) {
      T8DG_ASSERT (t8dg_face_dof_values_is_valid (mortar->fluxvalue_plus[isubface]));
      return mortar->fluxvalue_plus[isubface];
    }
  }
  T8DG_ABORT ("The element does not border the mortar");
}

void
t8dg_mortar_array_calculate_flux_dof1D (t8dg_mortar_array_t * mortar_array, t8dg_dof_values_t * dof_values,
                                        int icomp, t8dg_numerical_flux1D_fn numerical_flux, void *numerical_flux_data, double time)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface;
  t8dg_mortar_t      *mortar;
  double              ghost_exchange_time;

  ghost_exchange_time = -sc_MPI_Wtime ();
  t8dg_dof_values_ghost_exchange (dof_values);
  ghost_exchange_time += sc_MPI_Wtime ();
  mortar_array->ghost_exchange_time += ghost_exchange_time;

  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {

      /*TODO: num_faces eclass dependent */
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);

        /*for all faces check wether mortar is already allocated/computed */
        if (mortar == NULL) {
          mortar = t8dg_mortar_new (mortar_array->forest, itree, ielement, iface, mortar_array);
          t8dg_mortar_array_set_all_pointers (mortar_array, mortar);
        }
        if (!mortar->valid) {
          t8dg_mortar_calculate_flux_dof1D (mortar, dof_values, icomp, numerical_flux, numerical_flux_data, time);
        }
        T8DG_ASSERT (t8dg_mortar_array_get_mortar (mortar_array, idata, iface) != NULL);
      }
    }
  }
}

void
t8dg_mortar_array_calculate_linear_flux3D (t8dg_mortar_array_t * mortar_array, t8dg_dof_values_t * dof_values,
                                           t8dg_linear_flux3D_fn linear_flux, void *flux_data,
                                           t8dg_numerical_linear_flux3D_fn numerical_flux, void *numerical_flux_data, double time)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface;
  t8dg_mortar_t      *mortar;
  double              ghost_exchange_time;

  ghost_exchange_time = -sc_MPI_Wtime ();
  t8dg_dof_values_ghost_exchange (dof_values);
  ghost_exchange_time += sc_MPI_Wtime ();
  mortar_array->ghost_exchange_time += ghost_exchange_time;

  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);
  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);
    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {

      /*TODO: num_faces eclass dependent */
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);

        /*for all faces check wether mortar is already allocated/computed */
        if (mortar == NULL) {
          mortar = t8dg_mortar_new (mortar_array->forest, itree, ielement, iface, mortar_array);
          t8dg_mortar_array_set_all_pointers (mortar_array, mortar);
        }
        if (!mortar->valid) {
          t8dg_mortar_calculate_linear_flux3D (mortar, dof_values, linear_flux, flux_data, numerical_flux, numerical_flux_data, time);
        }
        T8DG_ASSERT (t8dg_mortar_array_get_mortar (mortar_array, idata, iface) != NULL);
      }
    }
  }
}

t8dg_mortar_array_t *
t8dg_mortar_array_new_empty (t8_forest_t forest, t8dg_local_values_t * local_values)
{
  int                 iface;
  t8_locidx_t         idata;
  t8dg_mortar_array_t *mortar_array = T8DG_ALLOC_ZERO (t8dg_mortar_array_t, 1);

  t8_forest_ref (forest);
  mortar_array->forest = forest;
  mortar_array->num_local_elements = t8_forest_get_num_element (forest);
  mortar_array->num_total_elements = mortar_array->num_local_elements + t8_forest_get_num_ghosts (forest);
  mortar_array->local_values = local_values;
  mortar_array->max_num_faces = t8dg_local_values_get_max_num_faces (local_values);

  /*TODO: split in max_num_faces and tree dependent actual num_faces */
  for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
    mortar_array->mortars[iface] = sc_array_new_count (sizeof (t8dg_mortar_t *), mortar_array->num_total_elements);

    for (idata = 0; idata < mortar_array->num_total_elements; idata++) {
      t8dg_mortar_array_set_mortar (mortar_array, idata, iface, NULL);
    }
  }
  return mortar_array;
}

void
t8dg_mortar_array_invalidate_all (t8dg_mortar_array_t * mortar_array)
{
  t8_locidx_t         itree, ielement, idata;
  t8_locidx_t         num_trees, num_elems_in_tree;
  int                 iface;
  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);

  t8dg_mortar_t      *mortar;

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);
        if (mortar != NULL) {
          mortar->valid = 0;
        }
      }
    }
  }
}

void
t8dg_mortar_array_destroy (t8dg_mortar_array_t ** pmortar_array)
{

  t8_locidx_t         itree, ielement, idata, iface;
  t8_locidx_t         num_trees, num_elems_in_tree;

  t8dg_mortar_t      *mortar;
  t8dg_mortar_array_t *mortar_array;
  mortar_array = *pmortar_array;

  num_trees = t8_forest_get_num_local_trees (mortar_array->forest);

  for (itree = 0, idata = 0; itree < num_trees; itree++) {

    num_elems_in_tree = t8_forest_get_tree_num_elements (mortar_array->forest, itree);

    for (ielement = 0; ielement < num_elems_in_tree; ielement++, idata++) {
      for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
        mortar = t8dg_mortar_array_get_mortar (mortar_array, idata, iface);
        if (mortar != NULL) {
          t8dg_mortar_array_set_all_pointers_to_NULL (mortar_array, mortar);
          t8dg_mortar_destroy (&mortar);
        }
      }

    }
  }
  for (iface = 0; iface < mortar_array->max_num_faces; iface++) {
    sc_array_destroy (mortar_array->mortars[iface]);
    mortar_array->mortars[iface] = 0;
  }
  t8_forest_unref (&mortar_array->forest);
  T8DG_FREE (mortar_array);
  *pmortar_array = NULL;
}

void
t8dg_mortar_array_apply_element_boundary_integral (t8dg_mortar_array_t * mortar_array,
                                                   t8_locidx_t itree, t8_locidx_t ielement, sc_array_t * element_result_dof)
{
  int                 iface, num_faces;

  sc_array_t         *face_flux_dof;
  sc_array_t         *face_dof;
  sc_array_t         *summand;

  t8_locidx_t         idata = t8dg_itree_ielement_to_idata (mortar_array->forest, itree, ielement);

  t8dg_functionbasis_t *functionbasis;
  functionbasis = t8dg_global_values_get_functionbasis (t8dg_local_values_get_global_values (mortar_array->local_values, itree, ielement));

  num_faces = t8dg_functionbasis_get_num_face_functionbasis (functionbasis);

  summand = t8dg_element_dof_values_duplicate (element_result_dof);
  t8dg_element_dof_values_set_zero (element_result_dof);

  for (iface = 0; iface < num_faces; iface++) {
    face_flux_dof = t8dg_mortar_array_get_oriented_flux (mortar_array, idata, iface);
    face_dof = t8dg_face_dof_values_duplicate (face_flux_dof);

    t8dg_local_values_apply_face_mass_matrix (mortar_array->local_values, itree, ielement, iface, face_flux_dof, face_dof);

    t8dg_functionbasis_transform_face_dof_to_element_dof (functionbasis, iface, face_dof, summand);

    T8DG_ASSERT (t8dg_element_dof_values_is_valid (summand));
    t8dg_element_dof_values_axpy (1, summand, element_result_dof);
    sc_array_destroy (face_dof);
  }
  sc_array_destroy (summand);
}

double
t8dg_mortar_array_get_ghost_exchange_time (t8dg_mortar_array_t * mortar_array)
{
  return mortar_array->ghost_exchange_time;
}
