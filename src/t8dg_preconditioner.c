#include <t8_forest.h>
#include <t8_forest/t8_forest_iterate.h>
#include "t8dg_adapt.h"
#include "t8dg_preconditioner.h"

#if T8_WITH_PETSC

/* Block-Preconditioner routines for matrix-free operations */
/* currently their name imply it is only block-jacobi which is implemented, but block-gauss-seidel is also possible, it uses the exact same structure and members; renaming these routines still has to be done */
extern PetscErrorCode JacobiShellPCCreate (t8dg_block_preconditioner_ctx_t **, t8dg_linear_advection_diffusion_problem_t *,
                                           t8dg_dof_values_t **, PetscInt *, int);
extern PetscErrorCode JacobiShellPCSetUp (PC);
extern PetscErrorCode JacobiShellPCApply (PC, Vec, Vec);
extern PetscErrorCode JacobiShellPCDestroy (PC);
extern PetscErrorCode MatMult_MF_Jacobi_Preconditioner (Mat, Vec, Vec);

/* Matrix-free PETSc Matrix-Routines for solving a multigrid coarse level problem, restricting the fine level problem onto the coarse level and vice versa */
extern PetscErrorCode MatMult_MF_Coarse_LVL (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_Prolongation (Mat, Vec, Vec);
extern PetscErrorCode MatMult_MF_Restriction (Mat, Vec, Vec);

/* This function initializes a selected preconditioner */
void
t8dg_precon_initialize_preconditioner (PC * pc, int selector, t8dg_precon_general_preconditioner_t * general_precon, void *problem,
                                       t8dg_dof_values_t ** problem_dofs, t8dg_time_derivation_matrix_application_t time_derivative,
                                       Mat * A, PetscInt * vec_global_index)
{
  general_precon->preconditioner_setup_time = -sc_MPI_Wtime ();
  /* Fill the general preconditioner with the selected preconditioning routine */
  switch (selector) {
  case 0:
    /* No preconditioning is selected */
    t8dg_precon_init_without_preconditioning (pc);
    break;
  case 1:
    /* Block-Jacobi-preconditioner is selected */
    t8dg_precon_init_jacobi (problem, problem_dofs, A, pc, &general_precon->jacobi_preconditioner_ctx, vec_global_index, selector);
    break;
  case 2:
    /* Two-Level-Multigrid-preconditioner is selected */
    t8dg_precon_init_two_level_mg (problem, problem_dofs, time_derivative, A, pc, t8dg_adapt_multigrid_coarsen_finest_level,
                                   &general_precon->smoother, &general_precon->smoother_pc, &general_precon->Restriction,
                                   &general_precon->A_coarse, &(general_precon->coarse_solver), &general_precon->coarse_pc,
                                   &general_precon->Prolongation, &general_precon->coarse_lvl, &general_precon->res_prol_ctx,
                                   &general_precon->cmat_ctx, vec_global_index);
    break;
  case 3:
    /* Block-Gauss-Seidel-preconditioner is selected */
    t8dg_precon_init_jacobi (problem, problem_dofs, A, pc, &general_precon->jacobi_preconditioner_ctx, vec_global_index, selector);
    break;
  default:
    /* Default equals no preconditioning */
    t8dg_precon_init_without_preconditioning (pc);
    break;
  }
  general_precon->preconditioner_setup_time += sc_MPI_Wtime ();
}

/* This function desctroys a preconditioner and frees allocated memory which was used during the preconditioning routine */
void
t8dg_precon_destroy_preconditioner (PC * pc, int selector, t8dg_precon_general_preconditioner_t * general_precon)
{
  PetscErrorCode      ierr;

  switch (selector) {
  case 0:
    /* No preconditioning was selected, therefore, nothing needs to be cleaned up */
    break;
  case 1:
    /* Block-Jacobi-preconditioner is selected  */
    t8dg_precon_destroy_block_preconditioner (general_precon->jacobi_preconditioner_ctx);
    ierr = JacobiShellPCDestroy (*pc);
    break;
  case 2:
    /* Two-Level-Multigrid preconditioner is selected */
    t8dg_precon_destroy_two_level_mg (&general_precon->smoother, &general_precon->Restriction, &general_precon->A_coarse,
                                      &(general_precon->coarse_solver), &general_precon->Prolongation, &general_precon->coarse_lvl,
                                      &general_precon->cmat_ctx);
    break;
  case 3:
    /* Block-Gauss-Seidel-preconditioner is selected  */
    t8dg_precon_destroy_block_preconditioner (general_precon->jacobi_preconditioner_ctx);
    ierr = JacobiShellPCDestroy (*pc);
    break;
  default:
    /* No preconditioning was selected, therefore, nothing needs to be cleaned up */
    break;
  }
}

/* Initializes the without-preconditioning option */
void
t8dg_precon_init_without_preconditioning (PC * pc)
{
  PetscErrorCode      ierr;
  ierr = PCSetType (*pc, PCNONE);
  CHKERRQ (ierr);
}

/**********************************************************************************************/
/************************ Beginning of Jacobi-Preconditioner-Routines *************************/
/**********************************************************************************************/
/****************************************** Take 2 ********************************************/

/* Initializes a Block-Jacobi preconditioner */
void
t8dg_precon_init_jacobi (void *problem, t8dg_dof_values_t ** problem_dofs, Mat * A, PC * pc,
                         t8dg_block_preconditioner_ctx_t ** precon_jacobi_ctx, PetscInt * vec_global_index, int selector)
{
  PetscErrorCode      ierr;

  /* Create a PCSHELL which resembles a Block-Preconditioner */
  ierr = PCSetType (*pc, PCSHELL);
  CHKERRQ (ierr);

  /* Create an application context */
  JacobiShellPCCreate (precon_jacobi_ctx, problem, problem_dofs, vec_global_index, selector);

  /* Assign the block preconditioner context to the preconditioner */
  ierr = PCShellSetContext (*pc, *precon_jacobi_ctx);
  CHKERRQ (ierr);

  /* Set the function that applies the preconditioner */
  ierr = PCShellSetApply (*pc, JacobiShellPCApply);
  CHKERRQ (ierr);

  /* Set the Destroy Function */
  ierr = PCShellSetDestroy (*pc, JacobiShellPCDestroy);
  CHKERRQ (ierr);

  /* Assign an appropriate name */
  if (selector == 1) {
    ierr = PCShellSetName (*pc, "Block-Jacobi-Preconditioner");
    CHKERRQ (ierr);
  }
  else if (selector == 3) {
    ierr = PCShellSetName (*pc, "Block-Gauss-Seidel-Preconditioner");
    CHKERRQ (ierr);
  }
  /* Set up the preconditioner */
  JacobiShellPCSetUp (*pc);
}

/* Creates a Jacobi preconditioner with corresponding 'preconditioning context' */
PetscErrorCode
JacobiShellPCCreate (t8dg_block_preconditioner_ctx_t ** precon_ctx, t8dg_linear_advection_diffusion_problem_t * problem,
                     t8dg_dof_values_t ** pdof_array, PetscInt * indexing, int selector)
{
  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;
  int                 iter;

  /* Allocate a new preconditioning context */
  ierr = PetscNew (&ctx);
  CHKERRQ (ierr);

  /* Assign values to the preconditiong context */
  (ctx->jacobi_ctx).current_point_in_time = t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem));
  (ctx->jacobi_ctx).problem = (void *) problem;;
  (ctx->jacobi_ctx).problem_dofs = pdof_array;
  (ctx->jacobi_ctx).problem_dofs_derivation = t8dg_dof_values_duplicate (*pdof_array);
  /* This function is declared within 't8dg_advect_diff_problem_cxx' */
  (ctx->jacobi_ctx).jacobi_preconditioner_advect_diff_application = t8dg_advect_diff_problem_block_precon_time_derivative_variant;
  (ctx->jacobi_ctx).timestep = t8dg_timestepping_data_get_time_step (t8dg_advect_diff_problem_get_time_data (problem));
  (ctx->jacobi_ctx).current_implicit_coefficient = 1.0; /* in the case of the implicit euler method */

  /* process-local wise */
  (ctx->jacobi_ctx).num_local_dofs =
    (size_t) t8dg_dof_get_num_local_elements (*pdof_array) * t8dg_dof_get_max_num_element_dof (*pdof_array);

  /* Whether block-jacobi (== 1) or block-gauss-seidel (==3) has been choosen */
  (ctx->jacobi_ctx).selection = selector;

  /* Create a local indexing scheme */
  ierr = PetscMalloc1 ((ctx->jacobi_ctx).num_local_dofs, (&((ctx->jacobi_ctx).local_indexing)));
  CHKERRQ (ierr);
  for (iter = 0; iter < (ctx->jacobi_ctx).num_local_dofs; ++iter) {
    ((ctx->jacobi_ctx).local_indexing)[iter] = iter;
  }

  /* Global indexng scheme for the initial implicit system resulting from an implicit timestepping method */
  ctx->global_indexing = indexing;

  /* Assign the new preconditioning context */
  *precon_ctx = ctx;

  return 0;
}

/* Sets up a Block preconditioner - calculates the preconditioning matrix */
PetscErrorCode
JacobiShellPCSetUp (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;

  ierr = PCShellGetContext (pc, (void **) &ctx);
  CHKERRQ (ierr);

  /* Declare the linear system which has to be solved in order to mimic the application of the inverse preconditioning-matrix */
  /* Create a local vector which will hold the result of the application of the inverse preconditioning matrix */
  ierr = VecCreateSeq (PETSC_COMM_SELF, (ctx->jacobi_ctx).num_local_dofs, &ctx->u);
  CHKERRQ (ierr);

  /* Create a local vector which holds the not yet preconditioned vector; it is going to be the right hand side of the linear system */
  ierr = VecCreateSeq (PETSC_COMM_SELF, (ctx->jacobi_ctx).num_local_dofs, &ctx->f);
  CHKERRQ (ierr);

  /* Create alocal matrix which applies the preconditioner to a vector */
  ierr =
    MatCreateShell (PETSC_COMM_SELF, (ctx->jacobi_ctx).num_local_dofs, (ctx->jacobi_ctx).num_local_dofs, (ctx->jacobi_ctx).num_local_dofs,
                    (ctx->jacobi_ctx).num_local_dofs, &ctx->jacobi_ctx, &ctx->M_jacobi);
  CHKERRQ (ierr);
  /* Set the matrix application of a multiplication in which this matrix takes place in */
  ierr = MatShellSetOperation (ctx->M_jacobi, MATOP_MULT, (void (*)(void)) MatMult_MF_Jacobi_Preconditioner);
  CHKERRQ (ierr);

  /* Create a krylov subspace method as solver regarding the linear system */
  ierr = KSPCreate (PETSC_COMM_SELF, &ctx->jacobi_preconditioner);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (ctx->jacobi_preconditioner, ctx->M_jacobi, ctx->M_jacobi);
  CHKERRQ (ierr);
  ierr = KSPGetPC (ctx->jacobi_preconditioner, &ctx->pc_jacobi_preconditioner);
  CHKERRQ (ierr);
  /* Choose no preconditioer for this GMRES */
  t8dg_precon_init_without_preconditioning (&ctx->pc_jacobi_preconditioner);
  ierr = KSPSetType (ctx->jacobi_preconditioner, KSPGMRES);
  CHKERRQ (ierr);
  /* Set convergence tolerances of the preconditioning solver */
  ierr = KSPSetTolerances (ctx->jacobi_preconditioner, 1.0e-7, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Set the solver up ready to use */
  ierr = KSPSetUp (ctx->jacobi_preconditioner);
  CHKERRQ (ierr);

  return 0;
}

/* Applies the Jacobi preconditioner to a vector */
PetscErrorCode
JacobiShellPCApply (PC pc, Vec in, Vec out)
{

  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;

  /* Get the preconditioning context */
  ierr = PCShellGetContext (pc, (void **) &ctx);
  CHKERRQ (ierr);

  /* Take the 'in' vector as the right hand side, so that the GMRES application mimics the application of the inverse of the preconditioning matrix */

  /* It is allowed to copy from a parallel to a sequential vec */
  /* Receive the 'in' vector as the right hand side */
  ierr = VecCopy (in, ctx->f);
  CHKERRQ (ierr);

  /* Solve the system in order to obtain the the application of the inverse preconditioning matrix; the result is stored in u */
  ierr = KSPSolve (ctx->jacobi_preconditioner, ctx->f, ctx->u);
  CHKERRQ (ierr);

  /* View whether or not the GMRES did converge */
  ierr = KSPConvergedRateView (ctx->jacobi_preconditioner, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ (ierr);

  /* Copy the result to the 'out' vector */
  ierr = VecCopy (ctx->u, out);
  CHKERRQ (ierr);

  return 0;
}

/* Frees the resources used by the preconditioning methods */
PetscErrorCode
JacobiShellPCDestroy (PC pc)
{
  PetscErrorCode      ierr;
  t8dg_block_preconditioner_ctx_t *ctx;

  /* Get the preconditioning context */
  ierr = PCShellGetContext (pc, (void **) &ctx);
  CHKERRQ (ierr);

  /* Free the preconditioning context */
  ierr = PetscFree (ctx);
  CHKERRQ (ierr);

  return 0;
}

/* Destroys/frees the allocated memory of the block preconditioner */
void
t8dg_precon_destroy_block_preconditioner (t8dg_block_preconditioner_ctx_t * ctx)
{
  PetscErrorCode      ierr;

  /* Free dof values */
  t8dg_dof_values_destroy (&((ctx->jacobi_ctx).problem_dofs_derivation));

  /* Free the allocated space of PETSc objects */
  ierr = PetscFree ((ctx->jacobi_ctx).local_indexing);
  CHKERRQ (ierr);
  ierr = MatDestroy (&ctx->M_jacobi);
  CHKERRQ (ierr);
  ierr = VecDestroy (&ctx->f);
  CHKERRQ (ierr);
  ierr = VecDestroy (&ctx->u);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ctx->jacobi_preconditioner);
  CHKERRQ (ierr);

}

/* Matrix-free application which resembles the application of the function describing the block jacobi preconditioning process */
PetscErrorCode
MatMult_MF_Jacobi_Preconditioner (Mat A, Vec in, Vec out)
{

  PetscErrorCode      ierr;
  t8dg_precon_block_matrix_ctx_t *jacobi_ctx;

  ierr = MatShellGetContext (A, &jacobi_ctx);
  CHKERRQ (ierr);

  /* Copy the petsc vec entries into the dof_values of the problem in order to compute their derivation */
  t8dg_precon_write_vec_to_dof (&in, *(jacobi_ctx->problem_dofs), jacobi_ctx->num_local_dofs);

  /* Calculate the time derivation of the degrees of freedom in the jacobi-preconditioning variant */
  (jacobi_ctx->jacobi_preconditioner_advect_diff_application) (*(jacobi_ctx->problem_dofs), jacobi_ctx->problem_dofs_derivation,
                                                               jacobi_ctx->current_point_in_time, jacobi_ctx->problem,
                                                               jacobi_ctx->selection);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-(jacobi_ctx->timestep * jacobi_ctx->current_implicit_coefficient), jacobi_ctx->problem_dofs_derivation,
                        *(jacobi_ctx->problem_dofs));

  /* Copy the dof_values back to a petsc_vector */
  t8dg_precon_write_dof_to_vec (*(jacobi_ctx->problem_dofs), &out, jacobi_ctx->local_indexing, jacobi_ctx->num_local_dofs);

  return 0;
}

/**********************************************************************************************/
/*************************** End of Jacobi-Preconditioner-Routines ****************************/
/**********************************************************************************************/

/**********************************************************************************************/
/****************************** Beginning of Multigrid-Routines *******************************/
/**********************************************************************************************/

/* Initializes a 2-Level Multigrid preconditioner */
void
t8dg_precon_init_two_level_mg (void *problem, t8dg_dof_values_t ** problem_dofs, t8dg_time_derivation_matrix_application_t time_derivative,
                               Mat * A_fine, PC * mg_pc, t8dg_corase_lvl_adapt_func_t coarsening_func, KSP * smoother, PC * smoother_pc,
                               Mat * Restriction, Mat * A_coarse, KSP * coarse_solver, PC * coarse_pc, Mat * Prolongation,
                               t8dg_mg_coarse_lvl_t * coarse_lvl, t8dg_mg_interpolating_ctx_t * res_prol_ctx,
                               t8dg_coarse_matrix_ctx_t * cmat_ctx, PetscInt * vec_global_index)
{
  PetscErrorCode      ierr;

  t8dg_debugf ("Initialization of a 2-Level Multigrid Preconditioner\n");

  /* Call to set up the contexts needed by the multigrid routines and to construct the coarse level from the initial forest of the given problem */
  t8dg_mg_set_up_two_lvl_precon ((t8dg_linear_advection_diffusion_problem_t *) problem, problem_dofs, time_derivative,
                                 coarse_lvl, coarsening_func, res_prol_ctx, cmat_ctx, vec_global_index);

  /* Create the matrix-free application of the problem on the coarse level */
  t8dg_mg_create_coarse_lvl_matrix (A_coarse, coarse_lvl, cmat_ctx);

  /* Create the restriction matrix which interpolates a Vector from the fine level onto the coarse level */
  t8dg_mg_create_restriction_matrix (Restriction, res_prol_ctx);

  /* Create the prolongation matrix which interpolates a Vector from the coarse level onto the fine level */
  t8dg_mg_create_prolongation_matrix (Prolongation, res_prol_ctx);

  /* Set Coarse Solver, Smoother, Matrices and Operators */
  /* Set the Multigrid Preconditioner */
  ierr = PCSetType (*mg_pc, PCMG);
  CHKERRQ (ierr);
  /* Sett the number of levels; 2 means one fine and one coarse level */
  ierr = PCMGSetLevels (*mg_pc, 2, NULL);
  CHKERRQ (ierr);
  /* Choose a V-Cycle routine for the Preconditioner */
  ierr = PCMGSetType (*mg_pc, PC_MG_MULTIPLICATIVE);
  CHKERRQ (ierr);
  ierr = PCMGSetCycleType (*mg_pc, PC_MG_CYCLE_V);
  CHKERRQ (ierr);
  /* Assign the coarse level solver */
  ierr = PCMGGetCoarseSolve (*mg_pc, coarse_solver);
  CHKERRQ (ierr);
  /* Set the operators for the coarse level solver and assign the corresponding matrix-free application */
  ierr = KSPSetOperators (*coarse_solver, *A_coarse, *A_coarse);
  CHKERRQ (ierr);
  /* Choose default values for convergence tolerances on the coarse level */
  ierr = KSPSetTolerances (*coarse_solver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  CHKERRQ (ierr);
  /* Choose a GMRES routine as the coarse level solver */
  ierr = KSPSetType (*coarse_solver, KSPGMRES);
  CHKERRQ (ierr);
  /* Extract the preconditioner on the coarse level */
  ierr = KSPGetPC (*coarse_solver, coarse_pc);
  CHKERRQ (ierr);
  /* Set no preconditioning on the coarse level */
  t8dg_precon_init_without_preconditioning (coarse_pc);
  /* Assign the smoother of the multigrid preconditioner */
  ierr = PCMGGetSmoother (*mg_pc, 1, smoother);
  CHKERRQ (ierr);
  /* The smoother performs on the initial problem/fine level and therefore, uses the initial matrix-application of the problem */
  ierr = KSPSetOperators (*smoother, *A_fine, *A_fine);
  CHKERRQ (ierr);
  /* Choose a GMRES routine as the smooter on the fine level */
  ierr = KSPSetType (*smoother, KSPGMRES);
  CHKERRQ (ierr);
  /* Extract the preconditioner of the smoother */
  ierr = KSPGetPC (*smoother, smoother_pc);
  CHKERRQ (ierr);
  /* Set no preconditioning within the smoothing iterations */
  t8dg_precon_init_without_preconditioning (smoother_pc);
  /* Amount of pre- and post-smoothing steps */
  ierr = PCMGSetNumberSmooth (*mg_pc, 3);
  CHKERRQ (ierr);
  /* Set the Restriction operator from the fine level onto the coarse level */
  ierr = PCMGSetRestriction (*mg_pc, 1, *Restriction);
  CHKERRQ (ierr);
  /* Set the Prolongation operator from the coarse level onto the fine level */
  ierr = PCMGSetInterpolation (*mg_pc, 1, *Prolongation);
  CHKERRQ (ierr);
  /* Set residual calculation - my be default since only Mat's and Vec's are used */
  ierr = PCMGSetResidual (*mg_pc, 1, PCMGResidualDefault, *A_fine);
  CHKERRQ (ierr);
  /* Provide work space vectors */
  /* Normally, 3 vectors are needed per multigrid levels, but PETSc manages their allocation and use on it's own */

  t8dg_debugf ("Multigrid components have been initialized\n");
}

/* Frees the allocated data needed by the two level multigrid preconditioner */
void
t8dg_precon_destroy_two_level_mg (KSP * smoother, Mat * Restriction, Mat * A_coarse, KSP * coarse_solver, Mat * Prolongation,
                                  t8dg_mg_coarse_lvl_t * coarse_lvl, t8dg_coarse_matrix_ctx_t * cmat_ctx)
{
  PetscErrorCode      ierr;

  /* Un-reference the coarsened forest describing the coarse level within in the multigrid preconditioner */
  t8_forest_unref (&coarse_lvl->forest_coarsened);

  /* Destroy the adapat_data, consisiting of allocated space for restriction and prolongation operations */
  t8dg_values_destroy_adapt_data (coarse_lvl->dg_values, &coarse_lvl->tmp_mortar_coarse_lvl);

  /* Destroy the used degrees of freedom needed by the restriction and prolongation routines */
  t8dg_dof_values_destroy ((coarse_lvl->dof_values_adapt));
  t8dg_dof_values_destroy (&(cmat_ctx->problem_dofs_derivation));
  t8dg_dof_values_destroy (&(coarse_lvl->adapt_data->dof_values));;
  t8dg_dof_values_destroy (&(coarse_lvl->adapt_data->dof_values_adapt));

  /* Destroy the PETSc Matrices used within the multigrid preconditioner */
  ierr = MatDestroy (A_coarse);
  CHKERRQ (ierr);
  ierr = MatDestroy (Restriction);
  CHKERRQ (ierr);
  ierr = MatDestroy (Prolongation);
  CHKERRQ (ierr);
  /* Destroy the PETSc krylov subspace methods/contexts */
  ierr = KSPDestroy (coarse_solver);
  CHKERRQ (ierr);
  ierr = KSPDestroy (smoother);
  CHKERRQ (ierr);
  /* Free the coarse level indexing scheme */
  ierr = PetscFree (coarse_lvl->global_indexing);
  CHKERRQ (ierr);
}

/* A function that sets up the components needed in order to perform a two level multigrid V-Cycle */
void
t8dg_mg_set_up_two_lvl_precon (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                               t8dg_time_matrix_application time_derivative, t8dg_mg_coarse_lvl_t * coarse_lvl,
                               t8dg_corase_lvl_adapt_func_t coarsening_func, t8dg_mg_interpolating_ctx_t * res_prol_ctx,
                               t8dg_coarse_matrix_ctx_t * cmat_ctx, PetscInt * fine_forest_indexing)
{

  /* Build a coarse mesh which is used within the multigrid preconditioner */
  t8dg_mg_construct_coarse_lvl (problem, pdof_array, coarse_lvl, coarsening_func);

  /* Assign some data regarding the coarse matrix context */
  cmat_ctx->coarse_lvl = coarse_lvl;
  //cmat_ctx->problem_dofs = *(coarse_lvl->dof_values);
  cmat_ctx->current_point_in_time = t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem));
  cmat_ctx->user_data = (void *) problem;;
  cmat_ctx->problem_dofs_derivation = t8dg_dof_values_duplicate (*(coarse_lvl->dof_values_adapt));
  cmat_ctx->time_derivative_func = time_derivative;
  cmat_ctx->timestep = t8dg_timestepping_data_get_time_step (t8dg_advect_diff_problem_get_time_data (problem));
  cmat_ctx->current_coefficient = 1.0;  /* in the case of the implicit euler method */

  /* Fill the interpolating_context with information about the coarse and fine level */
  /* Number of fine grid process-local degrees of freedom */
  res_prol_ctx->num_local_dofs_fine_grid =
    (size_t) (t8dg_dof_get_num_local_elements (*(coarse_lvl->dof_values)) * t8dg_dof_get_max_num_element_dof (*(coarse_lvl->dof_values)));

  /* Indexing vector for petsc entries on the initial/fine mesh */
  res_prol_ctx->fine_lvl_global_indexing = fine_forest_indexing;

  /* Number of coarse grid process-local degrees of freedom */
  res_prol_ctx->num_local_dofs_coarse_grid = coarse_lvl->num_local_dofs;

  /* Save the coarse level forest inside the restriction_prolongation context */
  res_prol_ctx->coarse_lvl = coarse_lvl;;

}

/* A function that constructs the coarse level for the multigrid preconditioner */
void
t8dg_mg_construct_coarse_lvl (t8dg_linear_advection_diffusion_problem_t * problem, t8dg_dof_values_t ** pdof_array,
                              t8dg_mg_coarse_lvl_t * coarse_lvl_mesh, t8dg_corase_lvl_adapt_func_t coarsening_func)
{
  PetscErrorCode      ierr;
  t8_gloidx_t         global_offset_to_first_local_elem;
  size_t              iter;

  t8dg_debugf ("The construction method of the coarse level has been called\n");

  /* Assign data to coarse level mesh context */
  coarse_lvl_mesh->problem = problem;
  coarse_lvl_mesh->dof_values = pdof_array;
  coarse_lvl_mesh->dg_values = *(t8dg_advect_diff_problem_get_dg_values (problem));
  coarse_lvl_mesh->adapt_data = *(t8dg_advect_diff_problem_get_adapt_data (problem));
  coarse_lvl_mesh->tmp_mortar_coarse_lvl = NULL;
  coarse_lvl_mesh->problem_forest = *(t8dg_advect_diff_problem_get_forest (problem));
  coarse_lvl_mesh->dof_values_adapt = t8dg_advect_diff_problem_get_dof_values_adapt (problem);

  /* Sets the current time within the adapt_data */
  t8dg_adapt_data_set_time (coarse_lvl_mesh->adapt_data,
                            t8dg_timestepping_data_get_current_time (t8dg_advect_diff_problem_get_time_data (problem)));
  /* Clone the current values of the problem into the adapt_data->dof_values (these dofs will be adapted/restricted/prolongated) */
  coarse_lvl_mesh->adapt_data->dof_values = t8dg_dof_values_clone (*(coarse_lvl_mesh->dof_values));

  /* If source and sink terms are considered */
  if (coarse_lvl_mesh->adapt_data->source_sink_fn != NULL) {
    t8dg_adapt_data_interpolate_source_fn (coarse_lvl_mesh->adapt_data);
  }

  /* Keep the original forest of the problem */
  t8_forest_ref (coarse_lvl_mesh->problem_forest);

  /* Initialize a coarsened forest conatining the elements of the coarser muligrid mesh */
  t8_forest_init (&(coarse_lvl_mesh->forest_coarsened));

  /* Set the user-data pointer of the new forest */
  t8_forest_set_user_data (coarse_lvl_mesh->forest_coarsened, coarse_lvl_mesh->adapt_data);

  /* Set the adapt function coarsening the forest */
  t8_forest_set_adapt (coarse_lvl_mesh->forest_coarsened, coarse_lvl_mesh->problem_forest, coarsening_func, 0);

  /* Ghost values are needed for solving the system on the coarser mesh */
  t8_forest_set_ghost (coarse_lvl_mesh->forest_coarsened, 1, T8_GHOST_FACES);

  /* Commit the pre-set forest, so that it will be adapted */
  t8_forest_commit (coarse_lvl_mesh->forest_coarsened);

  /* If source and sink terms were considered in former calculations regarding the adapted degrees of freedom, they are going to be destroyed, because they are not needed anymore */
  if (coarse_lvl_mesh->adapt_data->source_sink_fn != NULL) {
    t8dg_dof_values_destroy (&(coarse_lvl_mesh->adapt_data->source_sink_dof));
  }
  /* Allocate space for new local_values and ghost values regarding the coarsened forest */
  t8dg_values_mg_allocate_adapt (coarse_lvl_mesh->dg_values, coarse_lvl_mesh->forest_coarsened);

  /* Create and allocate new degrees of freedom for the coarsened forest; they are needed by coarse level solver */
  *(coarse_lvl_mesh->dof_values_adapt) =
    t8dg_dof_values_new (coarse_lvl_mesh->forest_coarsened, t8dg_values_get_global_values_array (coarse_lvl_mesh->dg_values));

  /* Create and allocate space for the degrees of freedom which result from the adaption of th adapt_data->dof_values */
  coarse_lvl_mesh->adapt_data->dof_values_adapt =
    t8dg_dof_values_new (coarse_lvl_mesh->forest_coarsened, t8dg_values_get_global_values_array (coarse_lvl_mesh->dg_values));;

  /* Calculate the number of local degrees of freedom on the coarse level mesh */
  coarse_lvl_mesh->num_local_dofs =
    (size_t) (t8dg_dof_get_num_local_elements (coarse_lvl_mesh->adapt_data->dof_values_adapt) *
              t8dg_dof_get_max_num_element_dof (coarse_lvl_mesh->adapt_data->dof_values_adapt));

  /* Allocate space for an array which holds the global indices of the degrees of freedom */
  ierr = PetscMalloc1 (coarse_lvl_mesh->num_local_dofs, &coarse_lvl_mesh->global_indexing);

  /* Compute the offset of the first local element concerning the forest and the amount of degrees of freedom per element */
  global_offset_to_first_local_elem = t8_forest_get_first_local_element_id (coarse_lvl_mesh->forest_coarsened);
  if (global_offset_to_first_local_elem != 0) {
    global_offset_to_first_local_elem =
      global_offset_to_first_local_elem * t8dg_dof_get_max_num_element_dof (*(coarse_lvl_mesh->dof_values_adapt));
  }

  /* Fill the array of global indices */
  for (iter = 0; iter < coarse_lvl_mesh->num_local_dofs; ++iter) {
    (coarse_lvl_mesh->global_indexing)[iter] = iter + global_offset_to_first_local_elem;
  }
  t8dg_debugf ("End of coarse level construction mehtod\n");
}

/* Matrix-free PETSc Matrix routine which mimics the application of the coarse level system matrix on a PETSc Vector */
PetscErrorCode
MatMult_MF_Coarse_LVL (Mat A, Vec in, Vec out)
{

  PetscErrorCode      ierr;
  t8dg_coarse_matrix_ctx_t *c_appctx;

  /* Get the matrix-free application context */
  ierr = MatShellGetContext (A, &c_appctx);
  CHKERRQ (ierr);

  /* Copy the petsc vec entries into the dof_values of the problem in order to compute their derivation */
  t8dg_precon_write_vec_to_dof (&in, *(c_appctx->coarse_lvl->dof_values), c_appctx->coarse_lvl->num_local_dofs);

  /* Calculate the time derivation of the degrees of freedom */
  (c_appctx->time_derivative_func) (*(c_appctx->coarse_lvl->dof_values), c_appctx->problem_dofs_derivation, c_appctx->current_point_in_time,
                                    c_appctx->user_data);

  /* calculate the sum of the initially passed degrees of freedom ('in' Vector) and their (scaled) derivation */
  t8dg_dof_values_axpy (-c_appctx->timestep * c_appctx->current_coefficient, c_appctx->problem_dofs_derivation,
                        *(c_appctx->coarse_lvl->dof_values));

  /* Copy the dof_values back to a petsc_vector */
  t8dg_precon_write_dof_to_vec (*(c_appctx->coarse_lvl->dof_values), &out, c_appctx->coarse_lvl->global_indexing,
                                c_appctx->coarse_lvl->num_local_dofs);

  return 0;
}

/* Matrix-free PETSc Matrix routine which interpolates/restricts a PETSc Vector regarding the fine level onto the coarse level problem */
PetscErrorCode
MatMult_MF_Restriction (Mat A, Vec in, Vec out)
{
  PetscErrorCode      ierr;
  t8dg_mg_interpolating_ctx_t *appctx;
  t8dg_mg_coarse_lvl_t *coarse_lvl;

  /* Get the mtarix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Get the coarse level mesh */
  coarse_lvl = appctx->coarse_lvl;

  /* Write the entries of the 'in' vector in adapt_data->dof_values, these are getting adapted/restricted */
  t8dg_precon_write_vec_to_dof (&in, coarse_lvl->adapt_data->dof_values, appctx->num_local_dofs_fine_grid);

  /* Calling iterate_replace executes the adaption/restriction */
  t8_forest_iterate_replace (coarse_lvl->forest_coarsened, coarse_lvl->problem_forest, t8dg_adapt_replace);

  /* Write the degrees of freedom to a petsc vec */
  t8dg_precon_write_dof_to_vec (coarse_lvl->adapt_data->dof_values_adapt, &out, coarse_lvl->global_indexing,
                                appctx->num_local_dofs_coarse_grid);

  /* Swap the initial problem properties to the coarse level problem */
  t8dg_dof_values_swap ((coarse_lvl->dof_values), (coarse_lvl->dof_values_adapt));
  t8dg_dof_values_swap (&(coarse_lvl->adapt_data->dof_values), &(coarse_lvl->adapt_data->dof_values_adapt));
  t8dg_values_mg_swap_instances_to_coarse_lvl (coarse_lvl->dg_values);

  return 0;
}

/* Matrix-free PETSc Matrix routine which interpolates/prolongates a PETSc Vector regarding the coarse level problem onto the fine level */
PetscErrorCode
MatMult_MF_Prolongation (Mat A, Vec in, Vec out)
{

  PetscErrorCode      ierr;
  t8dg_mg_interpolating_ctx_t *appctx;
  t8dg_mg_coarse_lvl_t *coarse_lvl;

  /* Get the mtarix-free application context */
  ierr = MatShellGetContext (A, &appctx);
  CHKERRQ (ierr);

  /* Get the coarse level mesh */
  coarse_lvl = appctx->coarse_lvl;

  /* Write the entries of the 'in' vector in adapt_data->dof_values, these are getting adapted/restricted */
  t8dg_precon_write_vec_to_dof (&in, coarse_lvl->adapt_data->dof_values, appctx->num_local_dofs_coarse_grid);

  /* Set the user_data of the forest; this is needed for the adaption */
  t8_forest_set_user_data (coarse_lvl->problem_forest, coarse_lvl->adapt_data);

  /* Allocate space for the prolongated degrees of freedom */
  t8dg_values_mg_allocate_adapt (coarse_lvl->dg_values, coarse_lvl->problem_forest);

  /* Calling iterate_replace executes the adaption/prolongation */
  t8_forest_iterate_replace (coarse_lvl->problem_forest, coarse_lvl->forest_coarsened, t8dg_adapt_replace);

  /* Write the degrees of freedom to a petsc vec */
  t8dg_precon_write_dof_to_vec (coarse_lvl->adapt_data->dof_values_adapt, &out, appctx->fine_lvl_global_indexing,
                                appctx->num_local_dofs_fine_grid);

  /* Swap the coarse level problem to the initial fine level state */
  t8dg_dof_values_swap ((coarse_lvl->dof_values), (coarse_lvl->dof_values_adapt));
  t8dg_dof_values_swap (&(coarse_lvl->adapt_data->dof_values), &(coarse_lvl->adapt_data->dof_values_adapt));
  t8dg_values_mg_swap_instances_to_fine_lvl (coarse_lvl->dg_values);

  /*Allocate for the next restriction step */
  t8dg_values_mg_allocate_adapt (coarse_lvl->dg_values, coarse_lvl->forest_coarsened);

  return 0;
}

/* This functions creates a matrix-free coarse level application needed by the GMRES solver on the coarse level within the multigrid preconditioner */
void
t8dg_mg_create_coarse_lvl_matrix (Mat * A_coarse, t8dg_mg_coarse_lvl_t * coarse_lvl, t8dg_coarse_matrix_ctx_t * cmat_ctx)
{
  PetscErrorCode      ierr;
  /* Create the restricted/coarsened matrix which used to solve the problem on the coarse mesh */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, coarse_lvl->num_local_dofs, coarse_lvl->num_local_dofs, PETSC_DETERMINE, PETSC_DETERMINE, cmat_ctx,
                    A_coarse);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) * A_coarse, "Coarse-Level-Matrix Application");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation on the coarse level mesh */
  ierr = MatShellSetOperation (*A_coarse, MATOP_MULT, (void (*)(void)) MatMult_MF_Coarse_LVL);
  CHKERRQ (ierr);
}

/* This functions creates a matrix-free routine to restrict the problem of the fine level onto the coarse level within the multigrid preconditioner */
void
t8dg_mg_create_restriction_matrix (Mat * Restriction, t8dg_mg_interpolating_ctx_t * res_prol_ctx)
{
  PetscErrorCode      ierr;
  /* Define a Restriction matrix which interpolates a vector from the fine level onto the coarse level */
  /* During the restriction the initial problem gets changed in order to be solved on the coarse level */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, res_prol_ctx->num_local_dofs_coarse_grid, res_prol_ctx->num_local_dofs_fine_grid, PETSC_DETERMINE,
                    PETSC_DETERMINE, res_prol_ctx, Restriction);
  CHKERRQ (ierr);
  ierr = PetscObjectSetName ((PetscObject) * Restriction, "Restriction Matrix - interpolating a vector from the fine to the coarse level");
  CHKERRQ (ierr);

  /* Define the (multiplicative) MatVec-Operation which resembles the application of the restriction matrix */
  ierr = MatShellSetOperation (*Restriction, MATOP_MULT, (void (*)(void)) MatMult_MF_Restriction);
  CHKERRQ (ierr);
}

/* This functions creates a matrix-free routine to prolongate the problem of the coarse level onto the fine level within the multigrid preconditioner */
void
t8dg_mg_create_prolongation_matrix (Mat * Prolongation, t8dg_mg_interpolating_ctx_t * res_prol_ctx)
{
  PetscErrorCode      ierr;
  /* Define a Prolongation matrix which interpolates a vector from the coarse level onto the fine level */
  /* During the prolongation the intital problem gets re-transformed from the coarse level properties to the initial/fine problem state */
  ierr =
    MatCreateShell (PETSC_COMM_WORLD, res_prol_ctx->num_local_dofs_fine_grid, res_prol_ctx->num_local_dofs_coarse_grid, PETSC_DETERMINE,
                    PETSC_DETERMINE, res_prol_ctx, Prolongation);
  CHKERRQ (ierr);
  ierr =
    PetscObjectSetName ((PetscObject) * Prolongation, "Prolongation Matrix - interpolating a vector from the coarse to the fine level");
  CHKERRQ (ierr);
  /* Define the (multiplicative) MatVec-Operation which resembles the application of the prolongation matrix */
  ierr = MatShellSetOperation (*Prolongation, MATOP_MULT, (void (*)(void)) MatMult_MF_Prolongation);
  CHKERRQ (ierr);
}

/**********************************************************************************************/
/********************************* End of Multigrid-Routines **********************************/
/**********************************************************************************************/

/* A function that writes a t8dg_dof_values_t to a PETSc Vector */
void
t8dg_precon_write_dof_to_vec (t8dg_dof_values_t * dofs, Vec * p_vec, PetscInt * indexing, size_t num_local_dofs)
{
  PetscErrorCode      ierr;
  /* Sets the values of a petsc vector to the entries of a t8dg_dof_values_t */
  ierr = VecSetValues (*p_vec, num_local_dofs, indexing, t8dg_dof_values_get_double_pointer (dofs, 0), INSERT_VALUES);
  CHKERRQ (ierr);
  /* Assemble the petsc vector */
  ierr = VecAssemblyBegin (*p_vec);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (*p_vec);
  CHKERRQ (ierr);
}

/* A function that writes a PETSc Vector to a t8dg_dof_values_t */
void
t8dg_precon_write_vec_to_dof (Vec * p_vec, t8dg_dof_values_t * dofs, size_t num_local_dofs)
{
  PetscErrorCode      ierr;
  PetscScalar        *vec_reader;
  double             *dof_pointer;
  int                 dof_iter;

  /* Retrieve the local part of a petsc vector */
  ierr = VecGetArrayRead (*p_vec, &vec_reader);
  CHKERRQ (ierr);

  /* Get a double pointer to the data inside the sc_array_t of t8dg_dof_values_t */
  dof_pointer = t8dg_dof_values_get_double_pointer (dofs, 0);
  /* Overwrite the dof values with the entries of a petsc vector */
  /* it is assumed (not checked!) that the local petsc vector is of length dof_values->num_local_elements */
  for (dof_iter = 0; dof_iter < num_local_dofs; ++dof_iter) {
    dof_pointer[dof_iter] = (double) vec_reader[dof_iter];
  }
  /* Restore the array entries in the petsc vector */
  ierr = VecRestoreArrayRead (*p_vec, &vec_reader);
  CHKERRQ (ierr);
}

/* Thus function updates the preconditioner within the DIRK methods due to the varying coefficients of the different stages */
void
t8dg_precon_dirk_update_preconditioner (int selector, t8dg_precon_general_preconditioner_t * preconditioner,
                                        double stage_related_a_coefficient, double stage_current_time)
{
  switch (selector) {
  case 0:
    /* No preconditioning was selected, therefore, no updating is needed */
    break;
  case 1:
    /* Block-Jacobi-Preconditioner is selected */
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_implicit_coefficient = stage_related_a_coefficient;
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_point_in_time = stage_current_time;
    break;
  case 2:
    /* Two-Level-Multigrid preconditioner is selected */
    (preconditioner->cmat_ctx).current_coefficient = stage_related_a_coefficient;
    (preconditioner->cmat_ctx).current_point_in_time = stage_current_time;
    break;
  case 3:
    /* Block-Gauss-Seidel-Preconditioner is selected */
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_implicit_coefficient = stage_related_a_coefficient;
    (preconditioner->jacobi_preconditioner_ctx->jacobi_ctx).current_point_in_time = stage_current_time;
    break;
  default:
    break;
  }
}

double
t8dg_precon_get_setup_time (t8dg_precon_general_preconditioner_t * preconditioner)
{
  return (preconditioner->preconditioner_setup_time);
}

#endif
