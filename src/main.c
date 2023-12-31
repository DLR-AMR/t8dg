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

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_vec.h>

#include "t8dg.h"
#include "t8dg_advect_diff_problem.h"
#include "t8dg_coarse_geometry.h"
#include "t8dg_flux.h"
#include "t8dg_flux_implementation.h"
#include "t8dg_common.h"
#include "t8dg_adapt.h"

#include <sc_options.h>
#include <sc.h>

static int
t8dg_check_options (int icmesh, int initial_cond_arg,
                    int uniform_level, int max_level, int min_level,
                    int number_LGL_points, double start_time, double end_time, double cfl, double delta_t, int time_steps, int time_order,
                    int vtk_freq, int adapt_freq, int adapt_arg, double diffusion_coefficient, int numerical_flux_arg, int source_sink_arg)
{
  if (!(icmesh >= 0 && icmesh <= 11))
    return 0;
  if (!(initial_cond_arg >= 0 && initial_cond_arg <= 17))
    return 0;
  if (!(uniform_level >= 0 && uniform_level <= 30))
    return 0;
  if (!(max_level >= uniform_level && max_level <= 30))
    return 0;
  if (!(min_level >= 0 && min_level <= uniform_level))
    return 0;
  if (!(number_LGL_points >= 1 && number_LGL_points <= MAX_LGL_NUMBER))
    return 0;
  if (!(start_time < end_time))
    return 0;
  if (!((cfl > 0 && cfl <= 1) || (cfl == 0 && (delta_t > 0 || time_steps > 0))))
    return 0;
  if (!(time_order >= 1 && time_order <= 4))
    return 0;
  if (!(vtk_freq >= 0))
    return 0;
  if (!(adapt_freq >= 0))
    return 0;
  if (!(adapt_arg >= 0 && adapt_arg <= 3))
    return 0;
  if (diffusion_coefficient < 0)
    return 0;
  if (!(numerical_flux_arg >= 0 && numerical_flux_arg <= 2))
    return 0;
  if (!(source_sink_arg >= 0 && source_sink_arg <= 1))
    return 0;
  return 1;
}

int
main (int argc, char *argv[])
{
  int                 mpiret;
  sc_options_t       *opt;
  char                help[BUFSIZ];
  int                 parsed, helpme;

  int                 initial_cond_arg;
  int                 uniform_level, max_level, min_level;
  int                 time_order;
  int                 number_LGL_points;
  int                 vtk_freq;
  int                 adapt_freq;
  int                 adapt_arg;
  int                 icmesh;
  int                 source_sink_arg;
  double              flow_velocity;
  double              cfl;
  double              start_time;
  double              end_time;
  double              delta_t;
  double              diffusion_coefficient;
  const char         *prefix;
  int                 numerical_flux_arg;
  int                 time_steps;
  int                 refine_error;
  /* brief help message */

  /* long help message */

  snprintf (help, BUFSIZ, "This program solves the advection diffusion equation. \n");
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
#ifdef T8_ENABLE_DEBUG
  t8dg_init (SC_LP_DEBUG);
  t8_init (SC_LP_DEBUG);
#else
  t8dg_init (SC_LP_ESSENTIAL);
  t8_init (SC_LP_ESSENTIAL);
#endif

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &helpme, "Display a short help message.");
  sc_options_add_int (opt, 'l', "uniform_level", &uniform_level, 3, "The uniform initial refinement level of the mesh. Default: 3");
  sc_options_add_int (opt, 'r', "max_level", &max_level, -1, "The maximum refinement level(>=0). Default: uniform_level");
  sc_options_add_int (opt, 'R', "min_level", &min_level, -1, "The number of refinement levels(>=0). Default: uniform_level");

  sc_options_add_int (opt, 'L', "LGL", &number_LGL_points, 2, "The number of LGL basis points/basisfunctions in 1D. Default: 2");
  sc_options_add_int (opt, 'o', "time_order", &time_order, 2,
                      "The order used for the runge Kutta timestepping (1<= order <=4). Default: 2");

  sc_options_add_int (opt, 'm', "cmesh", &icmesh, 0, "Choose cmesh. Default: 0\n" "\t\t0: line 1 tree\n" "\t\t1: line 3 trees\n"
                      "\t\t2: diagonal line more trees\n" "\t\t3: square\n" "\t\t4: square different size trees\n" "\t\t5: square moebius\n"
                      "\t\t6: moebius more tree\n" "\t\t7: parallelogram\n" "\t\t8: cube\n" "\t\t9: circle ring\n"
                      "\t\t10: square half periodic\n" "\t\t11: cylinder ring\n");
  sc_options_add_double (opt, 'c', "flow_velocity", &flow_velocity, 1.0, "The flow velocity. Default: 1.0");
  sc_options_add_double (opt, 'd', "diff_coeff", &diffusion_coefficient, 0, "The diffusion coefficient. Default: 0");

  sc_options_add_int (opt, 'i', "initial_cond", &initial_cond_arg, 0, "Choose initial condition function. Default: 0\n"
                      "\t\t0: constant function\n" "\t\t1: 1D hat function\n" "\t\t2: 1D step function\n"
                      "\t\t3: cosine product function(can be used in all dimensions)\n" "\t\t4: norm\n" "\t\t5: 2D hat\n"
                      "\t\t6: 2D circle step function\n" "\t\t7: 2D triangle step function\n" "\t\t8: 3D sphere step function\n"
                      "\t\t9: circle ring sphere step function\n" "\t\t10: circle ring sin angle\n" "\t\t11: cylinder ring sin product\n"
                      "\t\t12: cylinder ring smooth ball\n" "\t\t13: 1D smooth indicator\n" "\t\t14: 2D smooth indicator\n"
                      "\t\t15: 3D smooth indicator\n" "\t\t16: constant 0");

  sc_options_add_int (opt, 'a', "adapt_freq", &adapt_freq, 1, "The number of steps until adapt. Default: 1\t (0 means no adapt)");
  sc_options_add_int (opt, 'A', "adapt_fn", &adapt_arg, 0,
                      "Choose Adapt Function. Default: 0\n" "\t\t0: mass-crit\n" "\t\t1: relative minmax-crit\n"
                      "\t\t2: criterion for smoothed indicator functions\n"
                      "\t\t3: criterion for smoothed indicator functions on hypercubes\n");

  sc_options_add_int (opt, 'v', "vkt_freq", &vtk_freq, 1, "The number of steps until new vtk output. Default: 1\t (0 means no vtk)");
  sc_options_add_string (opt, 'p', "vtk_prefix", &prefix, "t8dg_advect_diff", "prefix to the vtk output");

  sc_options_add_double (opt, 'C', "CFL", &cfl, 1.0, "The CFL number used to determine the timestep. Default: 1.0");
  sc_options_add_double (opt, 't', "start_time", &start_time, 0.0, "The start time of the solve. Default: 0.0");
  sc_options_add_double (opt, 'T', "end_time", &end_time, 1.0, "The end time of the solve. Default: 1.0");
  sc_options_add_double (opt, 'D', "delta_t", &delta_t, 0, "If CFL=0 is given, you can enter the constant timestep. Default 0");
  sc_options_add_int (opt, 'S', "number of timesteps", &time_steps, 1, "Choose number of timesteps. Default: 1\n"
                      "\t\t0: no source sink\n" "\t\t1: 3D cylinder ring");

  sc_options_add_int (opt, 'n', "numerical_flux", &numerical_flux_arg, 0, "Choose numerical fluxes for diffusion:\n"
                      "\t\t0: central\n" "\t\t1: alternating");

  sc_options_add_int (opt, 's', "source_sink_fn", &source_sink_arg, 0, "Choose source/sink function. Default: 0\n"
                      "\t\t0: no source sink\n" "\t\t1: 3D cylinder ring");

  sc_options_add_switch (opt, 'E', "refine_error", &refine_error, "Refines the grid before calculating the errors");

  parsed = sc_options_parse (t8dg_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if (max_level == -1)
    max_level = uniform_level;
  if (min_level == -1)
    min_level = uniform_level;

  if (helpme) {
    /* display help message and usage */
    t8dg_global_essentialf ("%s\n", help);
    sc_options_print_usage (t8dg_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && t8dg_check_options (icmesh, initial_cond_arg, uniform_level, max_level, min_level, number_LGL_points,
                                              start_time, end_time, cfl, delta_t, time_steps, time_order, vtk_freq, adapt_freq, adapt_arg,
                                              diffusion_coefficient, numerical_flux_arg, source_sink_arg)) {
    t8dg_linear_advection_diffusion_problem_t *problem;
    problem =
      t8dg_advect_diff_problem_init_arguments (icmesh, uniform_level, number_LGL_points, initial_cond_arg, flow_velocity,
                                               diffusion_coefficient, start_time, end_time, cfl, delta_t, time_steps, time_order,
                                               min_level, max_level, adapt_arg, adapt_freq, prefix, vtk_freq,
                                               numerical_flux_arg, source_sink_arg, refine_error, sc_MPI_COMM_WORLD);

    t8dg_advect_diff_solve (problem);

    t8dg_advect_diff_problem_destroy (&problem);
  }
  else {
    /* wrong usage */
    t8dg_global_productionf ("\n\tERROR:Wrong usage.\n\n");
    sc_options_print_usage (t8dg_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
