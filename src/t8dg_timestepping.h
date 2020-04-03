/*
 * timestepping.h
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

/** @file t8dg_timestepping.h */
#ifndef SRC_TIMESTEPPING_H_
#define SRC_TIMESTEPPING_H_
#include <sc_containers.h>
#include "t8dg.h"

typedef struct t8dg_timestepping_data t8dg_timestepping_data_t;

typedef void        (*t8dg_time_matrix_application) (const sc_array_t * src, sc_array_t * dest, const double t,
                                                     const void *application_data);

/** implements a single step of runge kutta with a-values in the butcher-tableau only on the first minor diagonal
 * The time derivative application is given by time_derivative
 */
void
 
 
 
 
 
 
 
 t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                     t8dg_timestepping_data_t * time_data, sc_array_t ** pdof_array, void *user_data);

void                t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c);

t8dg_timestepping_data_t *t8dg_timestepping_data_new (int time_order, double start_time, double end_time, double cfl);

void                t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data);

double              t8dg_timestepping_data_get_current_time (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time);

double              t8dg_timestepping_data_get_end_time (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_cfl (t8dg_timestepping_data_t * time_data);

int                 t8dg_timestepping_data_get_time_order (t8dg_timestepping_data_t * time_data);

double              t8dg_timestepping_data_get_time_step (t8dg_timestepping_data_t * time_data);

void                t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t);

int                 t8dg_timestepping_data_is_endtime_reached (t8dg_timestepping_data_t * time_data);

#endif /* SRC_TIMESTEPPING_H_ */
