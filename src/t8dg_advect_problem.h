/*
 * t8dg_advect.h
 *
 *  Created on: Mar 22, 2020
 *      Author: lukas
 */
/** @file t8dg_advect_problem.h */

#ifndef SRC_T8DG_ADVECT_H_
#define SRC_T8DG_ADVECT_H_

#include <t8_cmesh.h>
#include <sc.h>

#include "t8dg_flux.h"
#include "t8dg_timestepping.h"
#include "t8dg_coarse_geometry.h"

typedef struct t8dg_linear_advection_problem t8dg_linear_advection_problem_t;

/* Enum for statistics. */
typedef enum
{
  ADVECT_ADAPT = 0,             /* adapt runtime */
  ADVECT_PARTITION,             /* partition runtime */
  ADVECT_PARTITION_DATA,        /* data partitioning runtime */
  ADVECT_BALANCE,               /* balance runtime */
  ADVECT_GHOST,                 /* ghost runtime */
  ADVECT_GHOST_EXCHANGE,        /* ghost exchange runtime */
  ADVECT_GHOST_WAIT,            /* ghost exchange waittime */
  ADVECT_REPLACE,               /* forest_iterate_replace runtime */
  ADVECT_IO,                    /* vtk runtime */
  ADVECT_INIT,                  /* initialization runtime */
  ADVECT_AMR,                   /* AMR runtime (adapt+partition+ghost+balance) including data exchange (partition/ghost) */
  ADVECT_SOLVE,                 /* solver runtime */
  ADVECT_TOTAL,                 /* overall runtime */

  ADVECT_PARTITION_PROCS,       /* number of processes sent to in partition */
  ADVECT_BALANCE_ROUNDS,        /* number of rounds in balance */
  ADVECT_GHOST_SENT,            /* number of ghosts sent to other processes */
  ADVECT_ELEM_AVG,              /* average global number of elements (per time step) */
  ADVECT_ERROR_INF,             /* l_infty error */
  ADVECT_ERROR_2,               /* L_2 error */
  ADVECT_VOL_LOSS,              /* The loss in volume in percent */

  ADVECT_NUM_STATS,             /* The number of statistics that we measure */
  ADVECT_NUM_TIME_STATS = ADVECT_TOTAL + 1      /* The number of time statistics that we only want to have counted once */
} advect_stats_t;

t8dg_linear_advection_problem_t *t8dg_advect_problem_init (t8_cmesh_t cmesh,
                                                           t8dg_coarse_geometry_t * coarse_geometry,
                                                           int dim,
                                                           t8dg_scalar_function_3d_time_fn u_initial,
                                                           t8dg_flux_t * flux,
                                                           int uniform_level, int max_level,
                                                           int number_LGL_points, t8dg_timestepping_data_t * time_data,
                                                           t8_forest_adapt_t adapt_fn, sc_MPI_Comm comm);

void                t8dg_advect_problem_init_elements (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_destroy (t8dg_linear_advection_problem_t ** pproblem);

int                 t8dg_advect_problem_endtime_reached (t8dg_linear_advection_problem_t * problem);

int                 t8dg_advect_problem_get_stepnumber (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_set_time_step (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_printdof (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_write_vtk (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_advance_timestep (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_adapt (t8dg_linear_advection_problem_t * problem, int measure_time);

void                t8dg_advect_problem_partition (t8dg_linear_advection_problem_t * problem, int measure_time);

void                t8dg_advect_problem_compute_and_print_stats (t8dg_linear_advection_problem_t * problem);

void                t8dg_advect_problem_accumulate_stat (t8dg_linear_advection_problem_t * problem, advect_stats_t stat, double value);

double              t8dg_advect_problem_l_infty_rel (const t8dg_linear_advection_problem_t * problem);

double              t8dg_advect_problem_l2_rel (const t8dg_linear_advection_problem_t * problem);

int                 t8dg_advect_gradient_adapt (t8_forest_t forest,
                                                t8_forest_t forest_from,
                                                t8_locidx_t which_tree,
                                                t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements,
                                                t8_element_t * elements[]);

int                 t8dg_advect_mass_adapt (t8_forest_t forest,
                                            t8_forest_t forest_from,
                                            t8_locidx_t which_tree,
                                            t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, int num_elements, t8_element_t * elements[]);

#endif /* SRC_T8DG_ADVECT_H_ */
