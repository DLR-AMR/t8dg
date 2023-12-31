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

#include <gtest/gtest.h>
#include "../src/t8dg.h"
#include "../src/t8dg_global_values.h"
#include "../src/t8dg_values.h"
#include "../src/t8dg_local_values.h"
#include "../src/t8dg_values.h"
#include "../src/t8dg_sc_array.h"
#include "../src/t8dg_coarse_geometry.h"
#include "../src/t8dg_common.h"
#include "../src/t8dg_cmesh.h"

#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <sc_mpi.h>

#ifdef SC_ENABLE_MPI
#include <mpi.h>
#endif

#include <tuple>



/* *INDENT-OFF* */
class ValuesChildInterpolationTest2D : public ::testing::Test
/* *INDENT-ON* */
{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 1);
    default_scheme = t8_scheme_new_default_cxx ();

    forest = t8_forest_new_uniform (cmesh, default_scheme, 0, 1, sc_MPI_COMM_WORLD);

    coarse_geometry = t8dg_coarse_geometry_new_2D_linear ();

    values = t8dg_values_new_LGL_hypercube (2, 2, coarse_geometry, forest);

  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_values_destroy (&values);
    t8_forest_unref (&forest);
  }

  t8dg_values_t      *values;
  t8dg_coarse_geometry_t *coarse_geometry;
  t8_forest_t         forest;
};

typedef struct t8dg_test_values_problem
{
  t8dg_values_t      *values;
  t8dg_dof_values_t  *dof_values;
  t8dg_dof_values_t  *dof_values_adapt;
} t8dg_test_values_problem_t;

static void
t8dg_test_values_replace_all_by_parents (t8_forest_t forest_old,
                                         t8_forest_t forest_new,
                                         t8_locidx_t itree,
                                         t8_eclass_scheme_c * ts,
                                         int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new)
{
  t8dg_test_values_problem_t *problem;
  int                 num_children;

  problem = (t8dg_test_values_problem_t *) t8_forest_get_user_data (forest_new);

  t8_element_t       *element = t8_forest_get_element_in_tree (forest_new, itree, first_ielem_new);
  num_children = ts->t8_element_num_children (element);

#ifndef T8DG_ENABLE_MPI
  ASSERT_EQ (num_children, num_elems_old);
  ASSERT_EQ (num_elems_new, 1);
#else
  ASSERT_EQ_MPI (num_children, num_elems_old);
  ASSERT_EQ_MPI (num_elems_new, 1);
#endif

  t8dg_values_set_element_adapt (problem->values, itree, first_ielem_new);
  t8dg_values_transform_child_dof_to_parent_dof (problem->values, problem->dof_values, problem->dof_values_adapt, itree,
                                                 num_children, first_ielem_old, first_ielem_new);

}

static void
t8dg_test_values_replace_all_by_children (t8_forest_t forest_old,
                                          t8_forest_t forest_new,
                                          t8_locidx_t itree,
                                          t8_eclass_scheme_c * ts,
                                          int num_elems_old, t8_locidx_t first_ielem_old, int num_elems_new, t8_locidx_t first_ielem_new)
{
  t8dg_test_values_problem_t *problem;
  int                 ichild, num_children;

  problem = (t8dg_test_values_problem_t *) t8_forest_get_user_data (forest_new);

  t8_element_t       *element = t8_forest_get_element_in_tree (forest_old, itree, first_ielem_old);
  num_children = ts->t8_element_num_children (element);

#ifndef T8DG_ENABLE_MPI
  ASSERT_EQ (num_children, num_elems_new);
  ASSERT_EQ (num_elems_old, 1);
#else
  ASSERT_EQ_MPI (num_children, num_elems_new);
  ASSERT_EQ_MPI (num_elems_old, 1);
#endif

  for (ichild = 0; ichild < num_children; ichild++) {
    t8dg_values_set_element_adapt (problem->values, itree, first_ielem_new + ichild);
    t8dg_values_transform_parent_dof_to_child_dof (problem->values, problem->dof_values, problem->dof_values_adapt, itree,
                                                   first_ielem_old, first_ielem_new + ichild, ichild);
  }
}

int
t8dg_test_values_refine_all (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t itree,
                             t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  return 1;
}

int
t8dg_test_values_coarsen_all (t8_forest_t forest,
                              t8_forest_t forest_from,
                              t8_locidx_t itree,
                              t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[])
{
  if (num_elements > 1) {
    return -1;
  }
  return 0;
}

TEST_F (ValuesChildInterpolationTest2D, lgl2_element)
{
  t8dg_global_values_t *global_values = t8dg_values_get_global_values_array (values)[T8_ECLASS_QUAD];

#ifndef T8DG_ENABLE_MPI
  ASSERT_EQ (t8dg_global_values_get_num_faces (global_values), 4);
#else
  ASSERT_EQ_MPI (t8dg_global_values_get_num_faces (global_values), 4);
#endif

  t8dg_dof_values_t  *dof_values_adapt;
  t8dg_dof_values_t  *dof_children_expected;
  t8dg_dof_values_t  *dof_values;

  double              parent_dof[4] = { 1, 2, 3, 5 };
  double              child_dof[16] = { 1, 1.5, 2, 2.75,
    1.5, 2, 2.75, 3.5,
    2, 2.75, 3, 4,
    2.75, 3.5, 4, 5
  };

  dof_values = t8dg_dof_values_new_data_local (forest, t8dg_values_get_global_values_array (values), parent_dof, 4);

  t8_forest_t         forest_adapt;

  t8dg_test_values_problem_t problem = { values, dof_values, NULL };

  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &problem);
  t8_forest_set_adapt (forest_adapt, forest, t8dg_test_values_refine_all, 0);
  t8_forest_commit (forest_adapt);

  t8dg_values_allocate_adapt (values, forest_adapt);

  dof_values_adapt = t8dg_dof_values_new (forest_adapt, t8dg_values_get_global_values_array (values));
  problem.dof_values_adapt = dof_values_adapt;
  t8_forest_iterate_replace (forest_adapt, forest, t8dg_test_values_replace_all_by_children);

  t8dg_values_cleanup_adapt (values);

  forest = forest_adapt;
  forest_adapt = NULL;

  dof_children_expected = t8dg_dof_values_new_data_local (forest, t8dg_values_get_global_values_array (values), child_dof, 16);

#ifndef T8DG_ENABLE_MPI
  ASSERT_TRUE (t8dg_dof_values_equal (dof_values_adapt, dof_children_expected));

#else
  ASSERT_TRUE_MPI (t8dg_dof_values_equal (dof_values_adapt, dof_children_expected));

#endif

  t8dg_dof_values_destroy (&dof_values);
  t8dg_dof_values_destroy (&dof_values_adapt);
  t8dg_dof_values_destroy (&dof_children_expected);
}

/*TODO: do a test for Face child interpolation*/
/*TOOD: do a test for projection*/


/* *INDENT-OFF* */
class PrecomputedValuesProjectionTest1D : public ::testing::TestWithParam<int>
/* *INDENT-ON* */
{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;
    int                 num_trees;
    sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &num_trees);
    cmesh = t8dg_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD, num_trees);
    default_scheme = t8_scheme_new_default_cxx ();

    forest = t8_forest_new_uniform (cmesh, default_scheme, 1, 1, sc_MPI_COMM_WORLD);    /*uniform level 1, gets coarsened */
    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

    int                 num_lgl = GetParam ();
    values = t8dg_values_new_LGL_hypercube (1, num_lgl, coarse_geometry, forest);

    dof_values = t8dg_dof_values_new (forest, t8dg_values_get_global_values_array (values));
  }
  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_values_destroy (&values);
    t8_forest_unref (&forest);
    t8dg_dof_values_destroy (&dof_values);
  }
  t8dg_values_t      *values;
  t8dg_coarse_geometry_t *coarse_geometry;
  t8_forest_t         forest, forest_adapt;
  t8dg_dof_values_t  *dof_values;
};

INSTANTIATE_TEST_CASE_P (lglRange, PrecomputedValuesProjectionTest1D, testing::Range (2, MAX_LGL_NUMBER + 1),);

TEST_P (PrecomputedValuesProjectionTest1D, const_one)
{
  int                 num_lgl, idof;
  t8_forest_t         forest_adapt;
  t8dg_dof_values_t  *dof_values_adapt;
  sc_array_t         *element_dof_values_parent;

  double             *tree_vertices;
  tree_vertices = t8_forest_get_tree_vertices (forest, 0);
  EXPECT_NE (tree_vertices, (double *) NULL);

  t8dg_test_values_problem_t problem = { values, dof_values, NULL };
  t8dg_dof_values_set_all_values (dof_values, 1);

  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data (forest_adapt, &problem);
  t8_forest_set_adapt (forest_adapt, forest, t8dg_test_values_coarsen_all, 0);
  t8_forest_commit (forest_adapt);

  t8dg_values_allocate_adapt (values, forest_adapt);

  dof_values_adapt = t8dg_dof_values_new (forest_adapt, t8dg_values_get_global_values_array (values));
  problem.dof_values_adapt = dof_values_adapt;
  t8_forest_iterate_replace (forest_adapt, forest, t8dg_test_values_replace_all_by_parents);

  t8dg_values_cleanup_adapt (values);

  forest = forest_adapt;
  forest_adapt = NULL;

#ifndef T8DG_ENABLE_MPI
  ASSERT_EQ (t8_forest_get_num_element (forest), 1);
#else
  ASSERT_EQ_MPI (t8_forest_get_num_element (forest), 1);
#endif
  element_dof_values_parent = t8dg_dof_values_new_element_dof_values_view (dof_values_adapt, 0, 0);
  num_lgl = GetParam ();

  for (idof = 0; idof < num_lgl; idof++) {
    EXPECT_NEAR (t8dg_element_dof_values_get_value (element_dof_values_parent, idof), 1, 1e-10);
  }
  t8dg_element_dof_values_destroy (&element_dof_values_parent);

  t8dg_dof_values_destroy (&dof_values_adapt);

}

/* *INDENT-OFF* */
class PrecomputedValuesL2norm1D : public ::testing::TestWithParam<std::tuple<int,std::tuple<t8dg_scalar_function_3d_time_fn,double>>>
/* *INDENT-ON* */

{
protected:
  /* *INDENT-OFF* */
  void SetUp () override
  /* *INDENT-ON* */
  {
    int                 num_lgl;
    int                 level = 4;      //dependent on parameter?
    t8_cmesh_t          cmesh;
    t8_scheme_cxx_t    *default_scheme;

    cmesh = t8_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD);
    default_scheme = t8_scheme_new_default_cxx ();
    forest = t8_forest_new_uniform (cmesh, default_scheme, level, 1, sc_MPI_COMM_WORLD);
    coarse_geometry = t8dg_coarse_geometry_new_1D_linear ();

    num_lgl = std::get < 0 > (GetParam ());

    values = t8dg_values_new_LGL_hypercube (1, num_lgl, coarse_geometry, forest);
    dof_values = t8dg_dof_values_new (forest, t8dg_values_get_global_values_array (values));
  }

  /* *INDENT-OFF* */
  void TearDown () override
  /* *INDENT-ON* */
  {
    t8dg_values_destroy (&values);
    t8dg_dof_values_destroy (&dof_values);
    t8_forest_unref (&forest);
  }

  t8dg_values_t      *values;
  t8dg_coarse_geometry_t *coarse_geometry;

  t8dg_dof_values_t  *dof_values;
  t8_forest_t         forest;
};

static double
t8dg_test_const_one (const double x[3], const double t, void *data)
{
  return 1;
}

static double
t8dg_test_sinx (const double x[3], const double t, void *data)
{
  return sin (2 * M_PI * x[0]);
}

static double
t8dg_test_expx (const double x[3], const double t, void *data)
{
  return exp (x[0]);
}

static const
  std::tuple <
t8dg_scalar_function_3d_time_fn, double >
function1 (t8dg_test_const_one, 1);
static const
  std::tuple <
t8dg_scalar_function_3d_time_fn, double >function2 (t8dg_test_sinx, sqrt (0.5));
static const
  std::tuple <
t8dg_scalar_function_3d_time_fn, double >function3 (t8dg_test_expx, sqrt ((exp (2) - 1) / 2));

/* *INDENT-OFF* */
INSTANTIATE_TEST_CASE_P (lgl_and_functions, PrecomputedValuesL2norm1D,
    ::testing::Combine (
	::testing::Range (2, 8),
	::testing::Values (function1, function2, function3)),);
/* *INDENT-ON* */

TEST_P (PrecomputedValuesL2norm1D, test_functions)
{
  double              time = 0;
  double              norm;

  t8dg_values_interpolate_scalar_function_3d_time (values, std::get < 0 > (std::get < 1 > (GetParam ())), time, NULL, dof_values);

  norm = t8dg_values_norm_l2 (values, dof_values, sc_MPI_COMM_WORLD);
  EXPECT_NEAR (norm, std::get < 1 > (std::get < 1 > (GetParam ())), exp (-4 * std::get < 0 > (GetParam ())));
}
