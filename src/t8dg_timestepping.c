/*
 * timestepping.c
 *
 *  Created on: Mar 14, 2020
 *      Author: lukas
 */

#include<sc_containers.h>
#include "t8dg_dof.h"
#include "t8dg.h"
#include "t8dg_timestepping.h"
#include "t8dg_timestepping_dirk_coefficients.h"
#include "t8dg_values.h"
#include "t8dg_advect_diff_problem.h"
#include "t8dg_preconditioner.h"

#define TRY_MG_PRECONDITIONER 0

struct t8dg_timestepping_data
{
  int                 time_order;/**< time order of the Runge-kutta timestepping*/
  double              delta_t;  /**< time step */
  double              t;        /**< current time */
  double              T;        /**< end time */
  double              cfl;      /**< cfl number*/
  int                 step_number;
  int                 use_implicit_timestepping;
  int                 preconditioner_selection;
};

#if T8_WITH_PETSC
/* Matrix-free routine that mimics the application of the matrix resulting from the Implicit/backwards Euler_Method */
//extern PetscErrorCode MatMult_MF_Impl_Euler (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_Impl_Euler (Mat, Vec, Vec);
#endif

#if T8_WITH_PETSC
/* Matrix-free routine that mimics the application of the matrices resultung (during the steps) of the second and third order DIRK methods */
extern PetscErrorCode MatMult_MF_DIRK (Mat, Vec, Vec);
#endif

#if T8_WITH_PETSC
/* Jacobi-Preconditioner for matrix-free operations */
extern PetscErrorCode JacobiShellPCCreate (t8dg_timestepping_precon_jacobi_ctx_t **, double, double);
extern PetscErrorCode JacobiShellPCSetUp (PC, Mat, Vec, t8dg_dof_values_t *, size_t);
extern PetscErrorCode JacobiShellPCApply (PC, Vec x, Vec y);
extern PetscErrorCode JacobiShellPCDestroy (PC);
#endif
/*pre computed butcher tableau values for rk with values only on first diagonal*/
double              rk1_b[1] = { 1 };

double              rk2_a[2] = { 1 };
double              rk2_b[2] = { 0.5, 0.5 };

double              rk3_a[3] = { 1. / 3, 2. / 3 };
double              rk3_b[3] = { 1. / 4, 0, 3. / 4 };

double              rk4_a[4] = { 0.5, 0.5, 1 };
double              rk4_b[4] = { 1. / 6, 1. / 3, 1. / 3, 1. / 6 };

void
t8dg_runge_kutta_fill_coefficients (int time_order, double **prk_a, double **prk_b, double **prk_c)
{
  switch (time_order) {
  case 1:
    *prk_a = NULL;
    *prk_b = rk1_b;
    break;
  case 2:
    *prk_a = rk2_a;
    *prk_b = rk2_b;
    break;
  case 3:
    *prk_a = rk3_a;
    *prk_b = rk3_b;
    break;
  case 4:
    *prk_a = rk4_a;
    *prk_b = rk4_b;
    break;
  default:
    break;
  }
  *prk_c = *prk_a;
}

void
t8dg_timestepping_runge_kutta_step (t8dg_time_matrix_application time_derivative,
                                    t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
  int                 istep;
  double             *rk_a, *rk_b, *rk_c;
  t8dg_dof_values_t  *dof_beginning;
  t8dg_dof_values_t  *dof_change;
  t8dg_dof_values_t  *dof_new;
  t8dg_dof_values_t  *dof_step;
  double              time_beginning, time_current, time_step;
  int                 time_order = t8dg_timestepping_data_get_time_order (time_data);

  time_beginning = t8dg_timestepping_data_get_current_time (time_data);
  time_current = time_beginning;
  time_step = t8dg_timestepping_data_get_time_step (time_data);

  t8dg_runge_kutta_fill_coefficients (time_order, &rk_a, &rk_b, &rk_c);

  dof_beginning = t8dg_dof_values_clone (*pdof_array);
  dof_new = t8dg_dof_values_duplicate (dof_beginning);
  dof_step = t8dg_dof_values_duplicate (dof_beginning);
  dof_change = t8dg_dof_values_duplicate (dof_beginning);

  time_derivative (*pdof_array, dof_change, time_current, user_data);
  t8dg_dof_values_axpyz (rk_b[0] * time_step, dof_change, dof_beginning, dof_new);

  for (istep = 0; istep < time_order - 1; istep++) {
    /*calculate the y-value for which the derivative needs to be evaluated
     * since a has only values on the first minor diagonal, only the k from the step before and the original y is needed*/

    t8dg_dof_values_axpyz (rk_a[istep] * time_step, dof_change, dof_beginning, dof_step);
    t8dg_dof_values_swap (pdof_array, &dof_step);
    /* calculate the derivative at the step time and y value */

    time_current = time_beginning + rk_c[istep] * time_step;
    t8dg_timestepping_data_set_current_time (time_data, time_current);
    time_derivative (*pdof_array, dof_change, time_current, user_data);
    /*add weighted summand to result */
    t8dg_dof_values_axpy (rk_b[istep + 1] * time_step, dof_change, dof_new);
  }

  t8dg_timestepping_data_set_current_time (time_data, time_beginning + time_step);
  t8dg_dof_values_swap (pdof_array, &dof_new);

  t8dg_dof_values_destroy (&dof_beginning);
  t8dg_dof_values_destroy (&dof_change);
  t8dg_dof_values_destroy (&dof_new);
  t8dg_dof_values_destroy (&dof_step);
}

t8dg_timestepping_data_t *
t8dg_timestepping_data_new_cfl (int time_order, double start_time, double end_time, double cfl, int use_implicit_timestepping,
                                int preconditioner_selection)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = cfl;
  time_data->step_number = 0;
  time_data->delta_t = -1;
  time_data->use_implicit_timestepping = use_implicit_timestepping;
  time_data->preconditioner_selection = preconditioner_selection;
  return time_data;
}

t8dg_timestepping_data_t *
t8dg_timestepping_data_new_constant_timestep (int time_order, double start_time, double end_time, double delta_t,
                                              int use_implicit_timestepping, int preconditioner_selection)
{
  T8DG_ASSERT (time_order > 0);
  t8dg_timestepping_data_t *time_data = T8DG_ALLOC (t8dg_timestepping_data_t, 1);
  time_data->time_order = time_order;
  time_data->t = start_time;
  time_data->T = end_time;
  time_data->cfl = 0;
  time_data->step_number = 0;
  time_data->delta_t = delta_t;
  time_data->use_implicit_timestepping = use_implicit_timestepping;
  time_data->preconditioner_selection = preconditioner_selection;
  return time_data;
}

void
t8dg_timestepping_data_destroy (t8dg_timestepping_data_t ** ptime_data)
{
  t8dg_timestepping_data_t *time_data = *ptime_data;
  time_data->time_order = -1;
  time_data->t = -1;
  time_data->T = -1;
  time_data->cfl = -1;
  time_data->delta_t = -1;
  time_data->use_implicit_timestepping = -1;
  time_data->preconditioner_selection = -1;
  T8_FREE (time_data);
  *ptime_data = NULL;
}

double
t8dg_timestepping_data_get_current_time (t8dg_timestepping_data_t * time_data)
{
  return time_data->t;
}

void
t8dg_timestepping_data_set_current_time (t8dg_timestepping_data_t * time_data, double current_time)
{
  time_data->t = current_time;
}

double
t8dg_timestepping_data_get_end_time (t8dg_timestepping_data_t * time_data)
{
  return time_data->T;
}

double
t8dg_timestepping_data_get_cfl (t8dg_timestepping_data_t * time_data)
{
  return time_data->cfl;
}

int
t8dg_timestepping_data_get_time_order (t8dg_timestepping_data_t * time_data)
{
  return time_data->time_order;
}

double
t8dg_timestepping_data_get_time_step (t8dg_timestepping_data_t * time_data)
{
  return time_data->delta_t;
}

void
t8dg_timestepping_data_set_time_step (t8dg_timestepping_data_t * time_data, double delta_t)
{
  if (delta_t < (time_data->T - time_data->t)) {
    time_data->delta_t = delta_t;
  }
  else {
    time_data->delta_t = time_data->T - time_data->t;
  }
}

int
t8dg_timestepping_data_is_endtime_reached (t8dg_timestepping_data_t * time_data)
{
  return !(time_data->t < time_data->T);
}

int
t8dg_timestepping_data_get_step_number (t8dg_timestepping_data_t * time_data)
{
  return time_data->step_number;
}

void
t8dg_timestepping_data_increase_step_number (t8dg_timestepping_data_t * time_data)
{
  time_data->step_number++;
}

double
t8dg_timestepping_data_get_time_left (t8dg_timestepping_data_t * time_data)
{
  return time_data->T - time_data->t;
}

/* Calculates the approximation of the next time step by using the implicit Euler Method in a matrix-free fashion */
void
t8dg_timestepping_implicit_euler (t8dg_time_matrix_application time_derivative,
                                  t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data,
                                  int preconditioner_selection)
{
#if T8_WITH_PETSC

  Mat                 A;
  Vec                 u, f;
  KSP                 ksp;
  PC                  pc;
  PetscErrorCode      ierr;
  t8dg_timestepping_impl_euler_ctx_t appctx;
  t8dg_precon_general_preconditioner_t preconditioner;

  /* Initialize the context needed by the matrix application */
  t8dg_timestepping_init_impl_euler_appctx (time_derivative, time_data, pdof_array, user_data, &appctx);

  /* Advance a time step */
  t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

         /********** Setting up and Solving the LS ****************/

  /* Create the right-hand-side of the system resulting from the implicit Euler-Method */
  ierr = VecCreateMPI (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, &f);
  CHKERRQ (ierr);

  /* Fill the vector with the current degrees of freedom resulting from the former time step */
  t8dg_precon_write_dof_to_vec (*(appctx.future_local_dofs), &f, appctx.global_indexing);

  /* Assign a proper name */
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side");
  CHKERRQ (ierr);

  /* Create an Approximation Vector storing the solution of the LS */
  /* Use the size and allocation similar to the right-hand-side (-> f) */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation");
  CHKERRQ (ierr);

  /* Setting up matrix-free Matrix */
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context needed by the matrix-free MatVec multiplication */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE,
                    (void *) &appctx, &A);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix of Imp Euler");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_Impl_Euler);
  CHKERRQ (ierr);

  /* Setting up KSP Solver */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);
  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);
  /* Select GMRES as Solver; default: Restarted GMRES(30) */
  /* If GMRES is used as smoother -> KSPFGMRES muste be used */
  ierr = KSPSetType (ksp, KSPFGMRES);
  CHKERRQ (ierr);
  /* Set error tolerances in the (outer) GMRES iteration */
  ierr = KSPSetTolerances (ksp, 1.0e-8, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  t8dg_debugf ("KSPFGMRES was selected\n");

  /* Initialize the preconditioner */
  t8dg_precon_initialize_preconditioner (&pc, preconditioner_selection, &preconditioner, user_data, appctx.future_local_dofs,
                                         time_derivative, &A, appctx.global_indexing);

  /* Use eventually further specified properties via the command line, otherwise it is just going to be set up */
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  t8dg_debugf ("\nFGMRES got called\n");

  /* Solve the Linear System */
  ierr = KSPSolve (ksp, f, u);
  CHKERRQ (ierr);

  t8dg_debugf ("\nFGMRES solve completed\n");

  /* View whether or not the GMRES did converge */
  ierr = KSPConvergedRateView (ksp, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

        /******************* The LS has been solved *********************/

  /* Write the solution from the PETSc Vector to degrees of freedom of the problem */
  t8dg_precon_write_vec_to_dof (&u, *pdof_array);

  /* Free used matrix-application related components */
  t8dg_timestepping_destroy_impl_euler_appctx (&appctx);

  /* Destroy the preconditioner */
  t8dg_precon_destroy_preconditioner (preconditioner_selection, &preconditioner);

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  //ierr = KSPDestroy (&ksp);
  //CHKERRQ (ierr);

  t8dg_debugf ("\nImplicit Euler-Method has been completed\n");
#else
  t8dg_global_productionf
    ("t8code/t8dg is currently not configured with PETSc.\n\nIn order to use this function t8dg needs to be configured with PETSc.\n");
#endif
}

#if T8_WITH_PETSC
/* Function that mimics the multiplication of the system matrix resulting from the implicit Euler_Method */
PetscErrorCode
MatMult_MF_Impl_Euler (Mat A, Vec in, Vec out)
{
  t8dg_debugf ("in matmultimpleuler\n");
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  PetscErrorCode      ierr;
  //PetscScalar        *current_local_approx;
  //int                 dof_iter;
  //double             *dof_values_ptr;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Write the entries into the degrees of freedom of the problem, so that they can be derived by the time_derivate function */
  t8dg_precon_write_vec_to_dof (&in, *(appctx->future_local_dofs));

  /* Calculate the time derivative (concerning the next time step) of the degrees of freedom passed to the MatMult (= 'in' Vector) and store their derivation inside future_local_dofs_derivative */
  (appctx->time_derivative_func) (*(appctx->future_local_dofs), appctx->future_local_dofs_derivative, appctx->next_time_point,
                                  appctx->user_data);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-(appctx->timestep), appctx->future_local_dofs_derivative, *(appctx->future_local_dofs));

  /* Write the result of the application of the system matrix of the implicit euler method into the out vector */
  t8dg_precon_write_dof_to_vec (*(appctx->future_local_dofs), &out, appctx->global_indexing);

  return 0;
}
#endif

/* Calculates the approximation of the next time step using a DIRK method (either DIRK(2,2) or DIRK(3,3) can be selected) */
void
t8dg_timestepping_dirk (t8dg_time_matrix_application time_derivative,
                        t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data, int num_order_stages,
                        int preconditioner_selection)
{
#if T8_WITH_PETSC
  T8DG_ASSERT ((num_order_stages == 2 || num_order_stages == 3));

  /* Declare the general DIRK coefficients */
  double              (*dirk_a_coeff)[num_order_stages][num_order_stages];
  double              (*dirk_b_coeff)[num_order_stages];
  double              (*dirk_c_coeff)[num_order_stages];

  int                 step_iter;
  int                 iter;
  PetscErrorCode      ierr;
  PetscInt            iteration_count[num_order_stages];
  t8dg_timestepping_dirk_ctx_t appctx;

  /* This array holds the pointer to solution calculated of the different steps */
  t8dg_dof_values_t  *local_dofs_step[num_order_stages];

  /* Initial point in time when the DIRK method got called */
  double              initial_time = t8dg_timestepping_data_get_current_time (time_data);

  /* Variables of the Linear system */
  Mat                 A;
  Vec                 u, f;
  KSP                 ksp;
  PC                  pc;

  /* Declaration Preconditioner */
  t8dg_precon_general_preconditioner_t preconditioner;

  /* Select whether the order 2 or order 3 DIRK method has been choosen */
  if (num_order_stages == 2) {
    dirk_a_coeff = &dirk22_a_coeff;
    dirk_b_coeff = &dirk22_b_coeff;
    dirk_c_coeff = &dirk22_c_coeff;
  }
  else if (num_order_stages == 3) {
    dirk_a_coeff = &dirk33_a_coeff;
    dirk_b_coeff = &dirk33_b_coeff;
    dirk_c_coeff = &dirk33_c_coeff;
  }

  /* Initialize the context needed by the matrix application of the DIRK Method */
  t8dg_timestepping_init_dirk_appctx (time_derivative, time_data, pdof_array, user_data, &appctx);

  /*********************** Setting up the Vectors and Matrix ***************************/

  /* Create the right-hand-side of the system resulting from the implicit Euler-Method */
  ierr = VecCreateMPI (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, &f);
  CHKERRQ (ierr);
  /* Fill the vector with the current degrees of freedom resulting from the former time step */
  t8dg_precon_write_dof_to_vec (*(appctx.local_derivation_degrees), &f, appctx.global_indexing);
  /* Assign aproper name */
  ierr = PetscObjectSetName ((PetscObject) f, "Right Hand Side");
  CHKERRQ (ierr);

  /* Create an Approximation Vector storing the solution of the LS */
  /* Use the size and allocation similar to the right-hand-side (-> f) */
  ierr = VecDuplicate (f, &u);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) u, "Approximation");
  CHKERRQ (ierr);

  /* Setting up the matrix-free Matrix */
  /* Create a matrix shell with local dimensions equal to the dimension of the Vec containing the process-local degrees of freedom and add an application context needed by the matrix-free MatVec multiplication */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, (PetscInt) appctx.num_local_dofs, (PetscInt) appctx.num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE,
                    (void *) &appctx, &A);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) A, "Matrix-Free Application-Matrix of DIRK methods");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation */
  ierr = MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) MatMult_MF_DIRK);
  CHKERRQ (ierr);

  /* Setting up KSP Solver */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  /* Set KSP: Set Operators; here: matrix that defines LS also serve as the matrix that definies the preconditioner */
  ierr = KSPSetOperators (ksp, A, A);
  CHKERRQ (ierr);
  /* Extract KSP and PC Context from  KSP Context */
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ (ierr);
  /* Select FGMRES as Solver; otherwise a multigrid preconditioner can not be applied */
  ierr = KSPSetType (ksp, KSPFGMRES);
  CHKERRQ (ierr);
  /* Set the convergence tolerances of the KSP solver */
  ierr = KSPSetTolerances (ksp, 1.0e-8, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Select and initialize a Preconditioner */
  t8dg_precon_initialize_preconditioner (&pc, preconditioner_selection, &preconditioner, user_data, appctx.local_derivation_degrees,
                                         time_derivative, &A, appctx.global_indexing);
  /* Use eventually further specified properties via the command line, otherwise it is just going to be set up */
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ (ierr);

  /************************* Creation of Vectors and Matrix finished *************************/
  /******************************** Solving the system begins ********************************/

  t8dg_debugf ("\nDIRK Method of order %d has been called\n", num_order_stages);

  /* In every RKV step a LS has to be solved in order to obtain an intermediate solution which is than further processed */
  for (step_iter = 0; step_iter < num_order_stages; ++step_iter) {

    /* Set the next time point of the upcoming RKV step */
    appctx.next_time_point = initial_time + (((*dirk_c_coeff)[step_iter]) * appctx.timestep);
    t8dg_timestepping_data_set_current_time (time_data, appctx.next_time_point);

    /* Set the current coefficient needed in the matrix-free application for the next step */
    appctx.dirk_current_a_coeff = (*dirk_a_coeff)[step_iter][step_iter];

    /* Fill the future_local_dofs_step with zeros, so that the axpy product will not be irritated */
    t8dg_dof_values_set_zero (appctx.future_local_dofs_step);

    /* Accumulate the (weighted) derivations of the degrees of freedom of the former stages */
    for (iter = 0; iter < step_iter; ++iter) {
      t8dg_dof_values_axpy ((appctx.timestep * ((*dirk_a_coeff)[step_iter][iter])), local_dofs_step[iter], appctx.future_local_dofs_step);
    }
    /* Add the (weighted) derivations of the degrees of freedom of the former stages to the right-hand side */
    if (step_iter > 0) {
      /* Update the right hand side of the linear system */
      t8dg_precon_write_dof_to_vec (appctx.future_local_dofs_step, &f, appctx.global_indexing);
    }

    /* Update the coarse level solver by assigning the current leading a_coeff */
    t8dg_precon_dirk_update_coarse_lvl_solver (preconditioner_selection, &preconditioner, appctx.dirk_current_a_coeff,
                                               appctx.next_time_point);

    /* Solve the Linear System */
    ierr = KSPSolve (ksp, f, u);
    CHKERRQ (ierr);

    /* Get the iteration count */
    ierr = KSPGetIterationNumber (ksp, &iteration_count[step_iter]);
    CHKERRQ (ierr);

    /* View whether or not the GMRES did converge */
    ierr = KSPConvergedRateView (ksp, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ (ierr);

    /* Write the solution of the stage into local_derivation_degrees */
    t8dg_precon_write_vec_to_dof (&u, *(appctx.local_derivation_degrees));

    /* Create a new empty t8dg_dof_values_t which will hold the derived solution of the stage */
    local_dofs_step[step_iter] = t8dg_dof_values_duplicate (*(appctx.local_derivation_degrees));

    /* Store the derivation of the degrees of freedom calculated in the current step in local_dofs_step */
    (appctx.time_derivative_func) (*(appctx.local_derivation_degrees), local_dofs_step[step_iter], appctx.next_time_point,
                                   appctx.user_data);

  }

  /*********************** The Linear system has been solved at all stages ***********************/

  /* Add the values of the stages to the former approximation to obtain the degrees of freedom of the next step */
  for (iter = 0; iter < num_order_stages; ++iter) {
    t8dg_dof_values_axpy ((appctx.timestep * (*dirk_b_coeff)[iter]), local_dofs_step[iter], appctx.current_local_dofs);
  }

  /* Swap the former degrees of freedom with the newly calculated degrees of freedom */
  t8dg_dof_values_swap (pdof_array, &(appctx.current_local_dofs));

  t8dg_debugf ("\nDIRK Method completed\n");

  /* Destroy the application context */
  t8dg_timestepping_destroy_dirk_appctx (&appctx);

  /* Destroy the preconditioner */
  t8dg_precon_destroy_preconditioner (preconditioner_selection, &preconditioner);

  /* Free the degrees of freedom obtained in the single stages */
  for (iter = 0; iter < num_order_stages; ++iter) {
    t8dg_dof_values_destroy (&local_dofs_step[iter]);
  }

  /* Free PETSc resources */
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);
  //ierr = KSPDestroy (&ksp);
  //CHKERRQ (ierr);

  if (num_order_stages == 2) {
    t8dg_debugf ("Iteration Count:\nStage 1: %d\nStage 2: %d\n", iteration_count[0], iteration_count[1]);
  }
  else if (num_order_stages == 3) {
    t8dg_debugf ("Iteration Count:\nStage 1: %d\nStage 2: %d\nStage 3: %d\n", iteration_count[0], iteration_count[1], iteration_count[2]);
  }
#else
  t8dg_global_productionf
    ("t8code/t8dg is currently not configured with PETSc.\n\nIn order to use this function t8dg needs to be configured with PETSc.\n");
#endif
}

#if T8_WITH_PETSC
/* Initializes the apllication context needed by the matrix-free application of the system matrix resulting from the implicit euler method */
void
t8dg_timestepping_init_impl_euler_appctx (t8dg_time_matrix_application time_derivative, t8dg_timestepping_data_t * time_data,
                                          t8dg_dof_values_t ** pdof_array, void *user_data, t8dg_timestepping_impl_euler_ctx_t * appctx)
{
  size_t              iter;
  PetscErrorCode      ierr;
  t8_gloidx_t         global_offset_to_first_local_elem;

  /* Store the current timestep */
  appctx->timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Store the next point in time which is going to be calculated */
  appctx->next_time_point = t8dg_timestepping_data_get_current_time (time_data) + (double) appctx->timestep;

  /* Store the time_derivative_function and the problem related data as well as the degrees of freedom */
  appctx->time_derivative_func = time_derivative;
  appctx->time_data = time_data;
  appctx->user_data = user_data;

  /* Store the degrees of freedom in the matrix-free application context (the current (original) degrees of freedom will be directly overwritten) */
  appctx->future_local_dofs = pdof_array;
  /* Is going to store the derivation of the degrees of freedom -> the evaluation of the time_derivative function */
  appctx->future_local_dofs_derivative = t8dg_dof_values_duplicate (*(appctx->future_local_dofs));

  /* Calculate the process-local number of degrees of freedom */
  appctx->num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (*(appctx->future_local_dofs)) *
              t8dg_dof_get_max_num_element_dof (*(appctx->future_local_dofs)));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (appctx->num_local_dofs, &appctx->global_indexing);

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (t8dg_dof_values_get_forest (*(appctx->future_local_dofs)));
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem = global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (*(appctx->future_local_dofs));
  }

  /* Fill the array of global indices */
  for (iter = 0; iter < appctx->num_local_dofs; ++iter) {
    (appctx->global_indexing)[iter] = iter + global_offset_to_first_local_elem;
  }

}
#endif

#if T8_WITH_PETSC
void
t8dg_timestepping_init_dirk_appctx (t8dg_time_matrix_application time_derivative, t8dg_timestepping_data_t * time_data,
                                    t8dg_dof_values_t ** pdof_array, void *user_data, t8dg_timestepping_dirk_ctx_t * appctx)
{
  int                 iter;
  PetscErrorCode      ierr;
  t8_gloidx_t         global_offset_to_first_local_elem;

  /* Store the time_derivative_function and the problem related data */
  appctx->time_derivative_func = time_derivative;
  appctx->time_data = time_data;
  appctx->user_data = user_data;
  appctx->timestep = t8dg_timestepping_data_get_time_step (time_data);

  /* Clone the current degrees of freedom */
  appctx->local_derivation_degrees = pdof_array;
  appctx->current_local_dofs = t8dg_dof_values_clone (*(appctx->local_derivation_degrees));

  /* Duplicate the other t8dg_dof_values_t */
  appctx->future_local_dofs_derivative = t8dg_dof_values_duplicate (*(appctx->local_derivation_degrees));
  appctx->future_local_dofs_step = t8dg_dof_values_duplicate (*(appctx->local_derivation_degrees));

  /* Calculate the process-local number of degrees of freedom */
  appctx->num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (*(appctx->local_derivation_degrees)) *
              t8dg_dof_get_max_num_element_dof (*(appctx->local_derivation_degrees)));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (appctx->num_local_dofs, &(appctx->global_indexing));

  /* compute the offset of the first local element concerning the forest and the amout of degrees of freedom per element */
  global_offset_to_first_local_elem =
    t8_forest_get_first_local_element_id (t8dg_dof_values_get_forest (*(appctx->local_derivation_degrees)));
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem =
      global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (*(appctx->local_derivation_degrees));
  }

  /* Fill the process-local array of global indices */
  for (iter = 0; iter < appctx->num_local_dofs; ++iter) {
    (appctx->global_indexing)[iter] = iter + global_offset_to_first_local_elem;
  }

}
#endif

#if T8_WITH_PETSC
/* Destroys the allocated space needed by the matrix application of the implicit euler method */
void
t8dg_timestepping_destroy_impl_euler_appctx (t8dg_timestepping_impl_euler_ctx_t * appctx)
{
  PetscErrorCode      ierr;

  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx->future_local_dofs_derivative));

  /* Free the PETSc indexing scheme */
  ierr = PetscFree (appctx->global_indexing);
  CHKERRQ (ierr);

}
#endif

#if T8_WITH_PETSC
void
t8dg_timestepping_destroy_dirk_appctx (t8dg_timestepping_dirk_ctx_t * appctx)
{
  PetscErrorCode      ierr;

  /* Free used resources */
  t8dg_dof_values_destroy (&(appctx->future_local_dofs_derivative));
  t8dg_dof_values_destroy (&(appctx->future_local_dofs_step));
  t8dg_dof_values_destroy (&(appctx->current_local_dofs));

  /* Free the PETSc indexing scheme */
  ierr = PetscFree (appctx->global_indexing);
  CHKERRQ (ierr);
}
#endif

#if T8_WITH_PETSC
PetscErrorCode
MatMult_MF_DIRK (Mat A, Vec in, Vec out)
{
  t8dg_timestepping_dirk_ctx_t *appctx;
  PetscErrorCode      ierr;
  PetscScalar        *current_local_approx;
  int                 dof_iter;
  double             *dof_values_ptr;

  /* Get the mtarix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Write the 'in' Vector to a t8dg_dof_values_t, so that the derivation of these coefficients can be calculated */
  t8dg_precon_write_vec_to_dof (&in, *(appctx->local_derivation_degrees));
  /* Get the derivation of the degrees of freedom of the current RKV step and store them in appctx->future_local_dofs_derivative */
  (appctx->time_derivative_func) (*(appctx->local_derivation_degrees), appctx->future_local_dofs_derivative, appctx->next_time_point,
                                  appctx->user_data);

  /* Subtract the derivation of the degrees from the dofs of the 'in' Vector */
  t8dg_dof_values_axpy (-(appctx->timestep * appctx->dirk_current_a_coeff), appctx->future_local_dofs_derivative,
                        *(appctx->local_derivation_degrees));

  /* Write the result of the matrix-application to the 'in' Vector in the 'out' Vector */
  t8dg_precon_write_dof_to_vec (*(appctx->local_derivation_degrees), &out, appctx->global_indexing);

  return 0;
}
#endif

/* Select either an implicit or an explicit Runge Kutta method */
void
t8dg_timestepping_choose_impl_expl_method (t8dg_time_matrix_application time_derivative,
                                           t8dg_timestepping_data_t * time_data, t8dg_dof_values_t ** pdof_array, void *user_data)
{
  if (time_data->use_implicit_timestepping == 0) {
    /* An explicit time stepping method has been chosen */
    t8dg_debugf ("Explicit RKV of order %d has been called.\n", time_data->time_order);
    t8dg_timestepping_runge_kutta_step (time_derivative, time_data, pdof_array, user_data);
  }
  else if (time_data->use_implicit_timestepping == 1) {
    /* An implicit time stepping method has been chosen */
    switch (time_data->time_order) {
    case 1:
      t8dg_debugf ("The implicit Euler method has been called.\n");
      t8dg_timestepping_implicit_euler (time_derivative, time_data, pdof_array, user_data, time_data->preconditioner_selection);
      break;
    case 2:
      t8dg_debugf ("The DIRK(2,2) method has been called.\n");
      t8dg_timestepping_dirk (time_derivative, time_data, pdof_array, user_data, 2, time_data->preconditioner_selection);
      break;
    case 3:
      t8dg_debugf ("The DIRK(3,3) method has been called.\n");
      t8dg_timestepping_dirk (time_derivative, time_data, pdof_array, user_data, 3, time_data->preconditioner_selection);
      break;
    default:
      t8dg_debugf ("An implicit DIRK method of order %d has been called, which is not implemented.\n", time_data->time_order);
    }
  }
}

/************** Preconditioners **************/
/* Jacobi Preconditioning routines */
#if T8_WITH_PETSC
/* Creates a Jacobi preconditioner with corresponding 'preconditioning context' */
PetscErrorCode
JacobiShellPCCreate (t8dg_timestepping_precon_jacobi_ctx_t ** pcshell, double timestep, double b_coeff)
{
  PetscErrorCode      ierr;
  t8dg_timestepping_precon_jacobi_ctx_t *pcctx;

  /* Create a new jacobi preconditioning context */
  ierr = PetscNew (&pcctx);
  CHKERRQ (ierr);
  /* Set Vec that represents the matrix diagonal */
  pcctx->matrix_diagonal = 0;
  pcctx->timestep = timestep;
  pcctx->current_b_coeff = b_coeff;
  /* return the preconditioning context */
  *pcshell = pcctx;

  return 0;
}

/* Sets up the Jacobi preconditioner - calculates main diagonal of the system matrix */
PetscErrorCode
JacobiShellPCSetUp (PC pc, Mat pmat, Vec u, t8dg_dof_values_t * current_dofs, size_t num_total_elements)
{
  PetscErrorCode      ierr;

  t8dg_timestepping_precon_jacobi_ctx_t *pcshell;
  Vec                 mat_diagonal;
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  t8dg_dof_values_t  *identity_vector;
  double             *dof_pointer_identity_vec;

  ierr = MatShellGetContext (pmat, &appctx);
  CHKERRQ (ierr);

  identity_vector = t8dg_dof_values_duplicate (current_dofs);

  ierr = VecDuplicate (u, &mat_diagonal);
  CHKERRQ (ierr);

  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  t8dg_linear_advection_diffusion_problem_t *problem = (t8dg_linear_advection_diffusion_problem_t *) appctx->user_data;

#if 1
  t8dg_advect_diff_problem_jacobi_precon (problem, current_dofs, identity_vector, pcshell->timestep, appctx->num_local_dofs,
                                          appctx->num_local_dofs);

  dof_pointer_identity_vec = t8dg_dof_get_double_pointer_to_array (identity_vector);

  ierr =
    VecSetValues (mat_diagonal, appctx->num_local_dofs, appctx->global_indexing, t8dg_dof_get_double_pointer_to_array (identity_vector),
                  INSERT_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (mat_diagonal);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (mat_diagonal);
  CHKERRQ (ierr);
  ierr = VecReciprocal (mat_diagonal);
  CHKERRQ (ierr);
#endif

  t8dg_debugf ("\nSetting up preconditioner finished\n");
  t8dg_debugf ("Preconditoner looks like\n*************\n");
  VecView (mat_diagonal, PETSC_VIEWER_STDOUT_WORLD);
  t8dg_debugf ("**************\n\n\n");

  //VecSet(mat_diagonal, 1.0);
  pcshell->matrix_diagonal = mat_diagonal;

  t8dg_dof_values_destroy (&identity_vector);
#if 0
#if 0
  //t8dg_timestepping_dirk_ctx_t *appctx;
  t8dg_timestepping_impl_euler_ctx_t *appctx;
  t8dg_dof_values_t  *dof_values_current_copy;
  t8dg_dof_values_t  *identity_vector;
  double             *dof_pointer;
  double             *dof_pointer_identity_vec;
#endif
  Vec                 mat_diagonal;
#if 0
  /* Get the matrix-free application context */
  ierr = MatShellGetContext (pmat, &appctx);
  CHKERRQ (ierr);
#endif
  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  /* Preallocate the mat_diagonal vector and an identity vector */
  ierr = VecDuplicate (u, &mat_diagonal);
  CHKERRQ (ierr);
  ierr = VecSet (mat_diagonal, 1.0);
  CHKERRQ (ierr);
#if 1
  dof_values_current_copy = t8dg_dof_values_clone (current_dofs);
  identity_vector = t8dg_dof_values_duplicate (current_dofs);
  dof_pointer = t8dg_dof_get_double_pointer_to_array (current_dofs);
  for (int i = 0; i < num_total_elements; ++i) {
    dof_pointer[i] = 0.0;
  }
  /* Calculate the diagonal of the application matrix */
  /******************* !!!!!! nur im seriellen Fall möglich !!!!!!!!! *********/
  dof_pointer_identity_vec = t8dg_dof_get_double_pointer_to_array (identity_vector);
  for (int i = 0; i < appctx->num_local_dofs; ++i) {
    dof_pointer[i] = 1.0;
    (appctx->time_derivative_func) (current_dofs, identity_vector, appctx->next_time_point, appctx->user_data);
    t8dg_debugf ("\n*****\nAbleitung wert von identity_vec auswertung an stelle %d ist %f\n*****\n", i, dof_pointer_identity_vec[i]);
    ierr =
      VecSetValue (mat_diagonal, i, 1.0 - ((pcshell->timestep * pcshell->current_b_coeff) * (dof_pointer_identity_vec[i])), INSERT_VALUES);
    CHKERRQ (ierr);
    ierr = VecAssemblyBegin (mat_diagonal);
    CHKERRQ (ierr);
    ierr = VecAssemblyEnd (mat_diagonal);
    CHKERRQ (ierr);
    dof_pointer[i] = 0.0;
  }
  t8dg_dof_values_swap (&current_dofs, &dof_values_current_copy);
  t8dg_dof_values_destroy (&current_dofs);
  t8dg_dof_values_destroy (&identity_vector);
  /* Calculate the reciprocal of the diagonal */
  ierr = VecReciprocal (mat_diagonal);
  CHKERRQ (ierr);
  //VecSet(mat_diagonal, 1.0);
  /* Store the preconditoner */
#endif
  pcshell->matrix_diagonal = mat_diagonal;
#if 0
  t8dg_debugf ("\nSetting up preconditioner finished\n");
  t8dg_debugf ("Preconditoner looks like\n*************\n");
  VecView (mat_diagonal, PETSC_VIEWER_STDOUT_WORLD);
  t8dg_debugf ("**************\n\n\n");
#endif
#endif
  return 0;
}

/* Applies the Jacobi preconditioner to a vector */
PetscErrorCode
JacobiShellPCApply (PC pc, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_timestepping_precon_jacobi_ctx_t *pcshell;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  /* Pointwise multiplication of the main diagonal to the 'out' vector, equals preconditioning with the diagonal matrix extracted from the system matrix */
  ierr = VecPointwiseMult (out, in, pcshell->matrix_diagonal);
  CHKERRQ (ierr);

  return 0;
}

/* Frees the resources used by the preconditioning methods */
PetscErrorCode
JacobiShellPCDestroy (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_timestepping_precon_jacobi_ctx_t *pcshell;

  /* Get the preconditoning context */
  ierr = PCShellGetContext (pc, (void **) &pcshell);
  CHKERRQ (ierr);

  /* Destroy the vector which stores the main diagonal entries */
  //ierr = VecDestroy(pcshell->matrix_diagonal); CHKERRQ(ierr);

  /* Free the prconditioning context */
  ierr = PetscFree (pcshell);
  CHKERRQ (ierr);

  return 0;
}

#endif

/*********************************************/
