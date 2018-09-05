/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// external
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

// acados
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_gnsf.h"

#include "acados/utils/external_function_generic.h"
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados/utils/timing.h"


#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// blasfeo
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"

// wt model
#include "examples/c/wt_model_nx6/wt_model.h"

// x0 and u for simulation
#include "examples/c/wt_model_nx6/u_x0.c"



int main()
{

	/************************************************
	* initialization
	************************************************/

	acados_timer test_timer;
	acados_tic(&test_timer);

	const int nx 		= 6;
	const int nu 		= 2;
	const int nz 		= 0;
	const int np 		= 1;
	const int nx1		= nx;
	const int nz1		= 0;
	const int ny 		= 5;  //nx + nu;
	const int nuhat 	= 0;  //nx + nu;
	const int n_out 	= 1;

	bool gnsf_init = true;

	int nsim = 1;

	int NF = nx + nu; // columns of forward seed

	double T = 1; // simulation time

	double *x_sim = malloc(sizeof(double)*nx*(nsim+1));

	double x_ref_sol[nx];
	double z_ref_sol[nz];
    double S_forw_ref_sol[nx*NF];
    double S_adj_ref_sol[NF];
    double S_alg_ref_sol[nz * (nx + nu)];

    double error[nx];
	double error_z[nz];
    double error_S_forw[nx*NF];
    double error_S_adj[NF];
	double error_S_alg[nz * (nx + nu)];

    double norm_error, norm_error_forw, norm_error_adj, norm_error_z, norm_error_sens_alg;
	double rel_error_x, rel_error_forw, rel_error_adj, rel_error_z, rel_error_alg;


	for (int ii = 0; ii < nx; ii++)
		x_sim[ii] = x_ref[ii];

	/************************************************
	* external functions (explicit model)
	************************************************/

	// expl_ode_fun
	external_function_param_casadi expl_ode_fun;
	expl_ode_fun.casadi_fun = &casadi_expl_ode_fun;
	expl_ode_fun.casadi_work = &casadi_expl_ode_fun_work;
	expl_ode_fun.casadi_sparsity_in = &casadi_expl_ode_fun_sparsity_in;
	expl_ode_fun.casadi_sparsity_out = &casadi_expl_ode_fun_sparsity_out;
	expl_ode_fun.casadi_n_in = &casadi_expl_ode_fun_n_in;
	expl_ode_fun.casadi_n_out = &casadi_expl_ode_fun_n_out;
	external_function_param_casadi_create(&expl_ode_fun, np);

	// expl_vde_for
	external_function_param_casadi expl_vde_for;
	expl_vde_for.casadi_fun = &casadi_expl_vde_for;
	expl_vde_for.casadi_work = &casadi_expl_vde_for_work;
	expl_vde_for.casadi_sparsity_in = &casadi_expl_vde_for_sparsity_in;
	expl_vde_for.casadi_sparsity_out = &casadi_expl_vde_for_sparsity_out;
	expl_vde_for.casadi_n_in = &casadi_expl_vde_for_n_in;
	expl_vde_for.casadi_n_out = &casadi_expl_vde_for_n_out;
	external_function_param_casadi_create(&expl_vde_for, np);

	// expl_vde_adj
	external_function_param_casadi expl_vde_adj;
	expl_vde_adj.casadi_fun = &casadi_expl_vde_adj;
	expl_vde_adj.casadi_work = &casadi_expl_vde_adj_work;
	expl_vde_adj.casadi_sparsity_in = &casadi_expl_vde_adj_sparsity_in;
	expl_vde_adj.casadi_sparsity_out = &casadi_expl_vde_adj_sparsity_out;
	expl_vde_adj.casadi_n_in = &casadi_expl_vde_adj_n_in;
	expl_vde_adj.casadi_n_out = &casadi_expl_vde_adj_n_out;
	external_function_param_casadi_create(&expl_vde_adj, np);

	/************************************************
	* external functions (implicit model)
	************************************************/

	// impl_ode_fun
	external_function_param_casadi impl_ode_fun;
	impl_ode_fun.casadi_fun = &casadi_impl_ode_fun;
	impl_ode_fun.casadi_work = &casadi_impl_ode_fun_work;
	impl_ode_fun.casadi_sparsity_in = &casadi_impl_ode_fun_sparsity_in;
	impl_ode_fun.casadi_sparsity_out = &casadi_impl_ode_fun_sparsity_out;
	impl_ode_fun.casadi_n_in = &casadi_impl_ode_fun_n_in;
	impl_ode_fun.casadi_n_out = &casadi_impl_ode_fun_n_out;
	external_function_param_casadi_create(&impl_ode_fun, np);

	// impl_ode_fun_jac_x_xdot
	external_function_param_casadi impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_fun = &casadi_impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_work = &casadi_impl_ode_fun_jac_x_xdot_work;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_sparsity_in;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_sparsity_out;
	impl_ode_fun_jac_x_xdot.casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_n_in;
	impl_ode_fun_jac_x_xdot.casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_n_out;
	external_function_param_casadi_create(&impl_ode_fun_jac_x_xdot, np);

	// impl_ode_jac_x_xdot_u
	external_function_param_casadi impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_fun = &casadi_impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_work = &casadi_impl_ode_jac_x_xdot_u_work;
	impl_ode_jac_x_xdot_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_sparsity_in;
	impl_ode_jac_x_xdot_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_sparsity_out;
	impl_ode_jac_x_xdot_u.casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_n_in;
	impl_ode_jac_x_xdot_u.casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_n_out;
	external_function_param_casadi_create(&impl_ode_jac_x_xdot_u, np);

	/************************************************
	* external functions (Generalized Nonlinear Static Feedback (GNSF) model)
	************************************************/
	// phi_fun
	external_function_param_casadi phi_fun;
	phi_fun.casadi_fun            = &casadi_phi_fun;
	phi_fun.casadi_work           = &casadi_phi_fun_work;
	phi_fun.casadi_sparsity_in    = &casadi_phi_fun_sparsity_in;
	phi_fun.casadi_sparsity_out   = &casadi_phi_fun_sparsity_out;
	phi_fun.casadi_n_in           = &casadi_phi_fun_n_in;
	phi_fun.casadi_n_out          = &casadi_phi_fun_n_out;
	external_function_param_casadi_create(&phi_fun, np);

	// phi_inc_dy
	external_function_param_casadi phi_fun_jac_y;
	phi_fun_jac_y.casadi_fun            = &casadi_phi_fun_jac_y;
	phi_fun_jac_y.casadi_work           = &casadi_phi_fun_jac_y_work;
	phi_fun_jac_y.casadi_sparsity_in    = &casadi_phi_fun_jac_y_sparsity_in;
	phi_fun_jac_y.casadi_sparsity_out   = &casadi_phi_fun_jac_y_sparsity_out;
	phi_fun_jac_y.casadi_n_in           = &casadi_phi_fun_jac_y_n_in;
	phi_fun_jac_y.casadi_n_out          = &casadi_phi_fun_jac_y_n_out;
	external_function_param_casadi_create(&phi_fun_jac_y, np);

	// phi_jac_y_uhat
	external_function_param_casadi phi_jac_y_uhat;
	phi_jac_y_uhat.casadi_fun                = &casadi_phi_jac_y_uhat;
	phi_jac_y_uhat.casadi_work               = &casadi_phi_jac_y_uhat_work;
	phi_jac_y_uhat.casadi_sparsity_in        = &casadi_phi_jac_y_uhat_sparsity_in;
	phi_jac_y_uhat.casadi_sparsity_out       = &casadi_phi_jac_y_uhat_sparsity_out;
	phi_jac_y_uhat.casadi_n_in               = &casadi_phi_jac_y_uhat_n_in;
	phi_jac_y_uhat.casadi_n_out              = &casadi_phi_jac_y_uhat_n_out;

	external_function_param_casadi_create(&phi_jac_y_uhat, np);

	// f_lo_fun_jac_x1k1uz
	external_function_param_casadi f_lo_fun_jac_x1k1uz;
	f_lo_fun_jac_x1k1uz.casadi_fun            = &casadi_f_lo_fun_jac_x1k1uz;
	f_lo_fun_jac_x1k1uz.casadi_work           = &casadi_f_lo_fun_jac_x1k1uz_work;
	f_lo_fun_jac_x1k1uz.casadi_sparsity_in    = &casadi_f_lo_fun_jac_x1k1uz_sparsity_in;
	f_lo_fun_jac_x1k1uz.casadi_sparsity_out   = &casadi_f_lo_fun_jac_x1k1uz_sparsity_out;
	f_lo_fun_jac_x1k1uz.casadi_n_in           = &casadi_f_lo_fun_jac_x1k1uz_n_in;
	f_lo_fun_jac_x1k1uz.casadi_n_out          = &casadi_f_lo_fun_jac_x1k1uz_n_out;
	external_function_param_casadi_create(&f_lo_fun_jac_x1k1uz, np);

	// get_matrices_fun
	external_function_casadi get_matrices_fun;
	get_matrices_fun.casadi_fun            = &casadi_get_matrices_fun;
	get_matrices_fun.casadi_work           = &casadi_get_matrices_fun_work;
	get_matrices_fun.casadi_sparsity_in    = &casadi_get_matrices_fun_sparsity_in;
	get_matrices_fun.casadi_sparsity_out   = &casadi_get_matrices_fun_sparsity_out;
	get_matrices_fun.casadi_n_in           = &casadi_get_matrices_fun_n_in;
	get_matrices_fun.casadi_n_out          = &casadi_get_matrices_fun_n_out;
	external_function_casadi_create(&get_matrices_fun);


/************************************************
* generate reference solution
************************************************/
	/* sim plan & config */

	// choose plan
	sim_solver_plan plan;

	plan.sim_solver = GNSF; // or IRK

	// create correct config based on plan
	sim_solver_config *config = sim_config_create(plan);

	/************************************************
	* sim dims
	************************************************/

	void *dims = sim_dims_create(config);
	config->set_nx(dims, nx);
	config->set_nu(dims, nu);
	config->set_nz(dims, nz);

	// GNSF -- set additional dimensions
	sim_gnsf_dims *gnsf_dim;
	if (plan.sim_solver == GNSF)
	{
		gnsf_dim = (sim_gnsf_dims *) dims;
		gnsf_dim->nx1   = nx1;
		gnsf_dim->nz1   = nz1;
		gnsf_dim->ny    = ny;
		gnsf_dim->nuhat = nuhat;
		gnsf_dim->n_out = n_out;
	}

	/************************************************
	* sim opts
	************************************************/

	sim_rk_opts *opts = sim_opts_create(config, dims);

	opts->ns = 10; // number of stages in rk integrator
	opts->num_steps = 1000; // number of integration steps
	opts->newton_iter = 5;
	opts->jac_reuse = false;
	opts->sens_adj = true;
	opts->sens_forw = true;


	/************************************************
	* sim in / out
	************************************************/

	sim_in *in = sim_in_create(config, dims);

	sim_out *out = sim_out_create(config, dims);

	in->T = T;

	/* set model */
	switch (plan.sim_solver)
	{
		case IRK:  // IRK
		{
			sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
			sim_set_model(config, in, "impl_ode_fun_jac_x_xdot",
					&impl_ode_fun_jac_x_xdot);
			sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
			break;
		}
		case GNSF:  // GNSF
		{
			// set model funtions
			sim_set_model(config, in, "phi_fun", &phi_fun);
			sim_set_model(config, in, "phi_fun_jac_y", &phi_fun_jac_y);
			sim_set_model(config, in, "phi_jac_y_uhat", &phi_jac_y_uhat);
			sim_set_model(config, in, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);

			// import model matrices
			external_function_generic *get_model_matrices =
					(external_function_generic *) &get_matrices_fun;
			gnsf_model *model = (gnsf_model *) in->model;
			sim_gnsf_import_matrices(gnsf_dim, model, get_model_matrices);
			break;
		}
		default :
		{
			printf("\nnot enough sim solvers implemented!\n");
			exit(1);
		}
	}

	// seeds forw
	for (int ii = 0; ii < nx * NF; ii++)
		in->S_forw[ii] = 0.0;
	for (int ii = 0; ii < nx; ii++)
		in->S_forw[ii * (nx + 1)] = 1.0;

	// seeds adj
	for (int ii = 0; ii < nx; ii++)
		in->S_adj[ii] = 1.0;
	for (int ii = 0; ii < nu; ii++)
		in->S_adj[ii+nx] = 0.0;

	/************************************************
	* sim solver
	************************************************/
	// print solver info
	// printf("\n ===  GENERATE REFERENCE SOLUTION USING IRK with %d stages and %d steps === \n", opts->ns, opts->num_steps);

	sim_solver *sim_solver = sim_create(config, dims, opts);
	int acados_return;

	if (plan.sim_solver == GNSF){  // for gnsf: perform precomputation
		gnsf_model *model = (gnsf_model *) in->model;
		sim_gnsf_precompute(config, gnsf_dim, model, opts,
					sim_solver->mem, sim_solver->work, in->T);
	}


	// to avoid unstable behavior introduce a small pi-controller for rotor speed tracking
	// double uctrl = 0.0;
	// double uctrlI = 0.0;
	// double kI = 1e-1;
	// double kP = 10;
	// double tmp, ctrlErr;



	for (int ii = 0; ii < 1; ii++) // upper bound was nsim, but for error computation we just compare first simulation
	{
		// update initial state
		for (int jj = 0; jj < nx; jj++)
			in->x[jj] = x_sim[ii*nx+jj];

		// // compute inputs
		// for (int jj = 0; jj < nu; jj++)
		// 	in->u[jj] = u_sim[ii*nu+jj];
		// tmp = in->u[1] - uctrl;
		// in->u[1] = tmp>0.0 ? tmp : 0.0;
		in->u[1] = 0.0;

		// update parameters
		switch (plan.sim_solver)
		{
			case IRK:  // IRK
			{
				impl_ode_fun.set_param(&impl_ode_fun, p_sim+ii*np);
				impl_ode_fun_jac_x_xdot.set_param(&impl_ode_fun_jac_x_xdot, p_sim+ii*np);
				impl_ode_jac_x_xdot_u.set_param(&impl_ode_jac_x_xdot_u, p_sim+ii*np);
				break;
			}
			case GNSF:  // GNSF
			{
				phi_fun.set_param(&phi_fun, p_sim+ii*np);
				phi_fun_jac_y.set_param(&phi_fun_jac_y, p_sim+ii*np);
				phi_jac_y_uhat.set_param(&phi_jac_y_uhat, p_sim+ii*np);
				f_lo_fun_jac_x1k1uz.set_param(&f_lo_fun_jac_x1k1uz, p_sim+ii*np);
				break;
			}
		}

		// d_print_mat(1, nx, in->x, 1);
		// d_print_mat(1, nu, in->u, 1);

		// execute simulation step with current input and state
		acados_return = sim_solve(sim_solver, in, out);
		if (acados_return != 0)
		{
			printf("error in sim solver\n");
			return ACADOS_FAILURE;
		}
		// extract state at next time step
		for (int jj = 0; jj < nx; jj++)
			x_sim[(ii+1)*nx+jj] = out->xn[jj];

		// update PI-controller
		// ctrlErr = x_ref[nx*(ii+1)] - x_sim[nx*(ii+1)];
		// uctrlI = uctrlI + kI*ctrlErr*T;
		// uctrl = kP*ctrlErr + uctrlI;

		// if (ii < nsim-1)
		// 	printf("\nii = %d, sim error = %e\n", ii, ctrlErr);
	}

	/************************************************
	* printing
	************************************************/
	printf("\nxn: \n");
	d_print_exp_mat(1, nx, &x_sim[nx], 1);  // print first result

	double *S_forw_out = NULL;
	if(opts->sens_forw){
		S_forw_out = out->S_forw;
		printf("\nS_forw_out: \n");
		d_print_exp_mat(nx, NF, S_forw_out, nx);
	}

    // store reference solution
    for (int jj = 0; jj < nx; jj++)
        x_ref_sol[jj] = out->xn[jj];

    for (int jj = 0; jj < nx*NF; jj++)
        S_forw_ref_sol[jj] = out->S_forw[jj];

    for (int jj = 0; jj < NF; jj++)
        S_adj_ref_sol[jj] = out->S_adj[jj];

    for (int jj = 0; jj < nz; jj++)
        z_ref_sol[jj] = out->zn[jj];

    for (int jj = 0; jj < nz*NF; jj++)
        S_alg_ref_sol[jj] = out->S_algebraic[jj];

	// compute one norms
    double norm_x_ref, norm_S_forw_ref, norm_S_adj_ref, norm_z_ref, norm_S_alg_ref = 0;

    norm_x_ref = onenorm(nx, 1, x_ref_sol);
    norm_S_forw_ref = onenorm(nx, nx + nu, S_forw_ref_sol);
    norm_S_adj_ref = onenorm(1, nx + nu, S_adj_ref_sol);
    norm_z_ref = onenorm(nz, 1, z_ref_sol);
    norm_S_alg_ref = onenorm(nz, nx + nu, S_alg_ref_sol);

	free(config);
    free(dims);
    free(opts);

    free(in);
    free(out);
    free(sim_solver);

/************************************************
* numerical experiment
************************************************/
	int n_executions = 10;

	bool jac_reuse 	= false;
	bool sens_forw 	= true;
	bool sens_adj  	= true;
	bool output_z  	= false;
	bool sens_alg  	= false;
	bool sens_hess  = false;

	int max_num_stages = 9;
	int min_num_stages = 1;
	int stages_in_experiment = max_num_stages - min_num_stages;

	int steps_in_experiment = 7;
	int steps_array[steps_in_experiment];
	steps_array[0] = 1;
	steps_array[1] = 2;
	steps_array[2] = 5;
	steps_array[3] = 10;
	steps_array[4] = 20;
	steps_array[5] = 50;
	steps_array[6] = 100;

	int min_newton = 1;
	int max_newton = 4;
	int newton_in_experiment = max_newton - min_newton;

	int num_experiments = steps_in_experiment * stages_in_experiment * newton_in_experiment;

	/* arrays for numerical experiment */
	// options
	double experiment_num_stages[num_experiments];
	double experiment_num_steps[num_experiments];
	double experiment_solver[num_experiments];
	double experiment_jac_reuse[num_experiments];
	double experiment_newton_it[num_experiments];
	double experiment_sens_forw[num_experiments];
	double experiment_sens_adj[num_experiments];
	double experiment_output_z[num_experiments];
	double experiment_sens_alg[num_experiments];
	// errors
	double experiment_error_sim[num_experiments];
	double experiment_error_sens_forw[num_experiments];
	double experiment_error_sens_adj[num_experiments];
	double experiment_error_z[num_experiments];
	double experiment_error_sens_alg[num_experiments];
	// timings
	double experiment_cpu_time[num_experiments];
	double experiment_ad_time[num_experiments];
	double experiment_lss_time[num_experiments];
	double experiment_la_time[num_experiments];

	// store current times
	double cpu_time[nsim];
	double lss_time[nsim];
	double ad_time[nsim];
	double la_time[nsim];

	double cpu_time_experiment;
	double lss_time_experiment;
	double la_time_experiment;
	double ad_time_experiment;

	int number_sim_solvers = 4;
	for (int nss = 2; nss < number_sim_solvers; nss++)
	{
		/************************************************
		* nss: number of sim_solver
			1 : ERK
			2 : IRK
			3 : GNSF
		************************************************/
		printf("\n ===  USING SOLVER NUMBER %d === \n",nss);
		/* initialize experiment times  */ 
		for (int ii = 0; ii < num_experiments; ii++) {
			// with 1e15 such that minimum can be taken easily
			experiment_cpu_time[ii] = 1e15;
			experiment_ad_time[ii] = 1e15;
			experiment_lss_time[ii] = 1e15;
			experiment_la_time[ii] = 1e15;		
		}
		for (int i_execution = 0; i_execution < n_executions; i_execution++) {
			printf("**************************************\n");
			printf("**************************************\n");
			printf("**** TEST EXECUTION NO   %d  *******  \n", i_execution);
			printf("**************************************\n");
			printf("**************************************\n");			
			for (int num_stages = min_num_stages; num_stages < max_num_stages; num_stages++) {
				
			// print experiment info
			printf("num_stages %d \n", num_stages);

				for (int tested_steps = 0; tested_steps < steps_in_experiment; tested_steps++) {
					for (int newton_iter = min_newton; newton_iter < max_newton; newton_iter++){
					/* sim plan & config */
						sim_solver_plan plan;
						switch (nss)
						{
							case 1:
								plan.sim_solver = ERK;
								break;

							case 2:
								plan.sim_solver = IRK;
								break;

							case 3:
								plan.sim_solver = GNSF;
								break;

							default :
								printf("\nnot enough sim solvers implemented!\n");
								exit(1);

						}

						// create correct config based on plan
						sim_solver_config *config = sim_config_create(plan);



					/* sim dims */
						void *dims = sim_dims_create(config);
						config->set_nx(dims, nx);
						config->set_nu(dims, nu);
						config->set_nz(dims, nz);

						// GNSF -- set additional dimensions
						sim_gnsf_dims *gnsf_dim;
						if (plan.sim_solver == GNSF)
						{
							gnsf_dim = (sim_gnsf_dims *) dims;
							gnsf_dim->nx1   = nx1;
							gnsf_dim->nz1   = nz1;
							gnsf_dim->ny    = ny;
							gnsf_dim->nuhat = nuhat;
							gnsf_dim->n_out = n_out;
						}


					/* sim options */

						void *opts_ = sim_opts_create(config, dims);
						sim_rk_opts *opts = (sim_rk_opts *) opts_;
						config->opts_initialize_default(config, dims, opts);

						opts->newton_iter = newton_iter;        // number of newton iterations per integration step
						opts->ns                = num_stages;   // number of stages in rk integrator
						int num_steps = steps_array[tested_steps];
						opts->num_steps         = num_steps;    // number of steps

						opts->jac_reuse = jac_reuse;        	// jacobian reuse
						opts->sens_forw         = sens_forw;
						opts->sens_adj          = sens_adj;
						opts->output_z          = output_z;
						opts->sens_algebraic    = sens_alg;
						opts->sens_hess    		= sens_hess;

					/* sim in / out */

						sim_in *in = sim_in_create(config, dims);
						sim_out *out = sim_out_create(config, dims);

						in->T = T;

					/* set model */
						switch (plan.sim_solver)
						{
							case ERK:
							{
								sim_set_model(config, in, "expl_ode_fun", &expl_ode_fun);
								sim_set_model(config, in, "expl_vde_for", &expl_vde_for);
								sim_set_model(config, in, "expl_vde_adj", &expl_vde_adj);
								break;
							}
							case IRK:  // IRK
							{
								sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
								sim_set_model(config, in, "impl_ode_fun_jac_x_xdot",
										&impl_ode_fun_jac_x_xdot);
								sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
								/* initialize integration variables to be equivalent to GNSF */
								// xdot1_z1_0 = (E\(A*x0(1:nx1)+ B * u0 + c))'

								// xdot1_z1_0 =
								//   -3.290336833904700e-03                         0     1.353969828015453e+00     3.228986350219792e-04    -2.067303680241406e+01    -2.005255049361735e+02
								if (gnsf_init){
									in->xdot[0] = -3.290336833904700e-03;
									in->xdot[1] = 0;
									in->xdot[2] = 1.353969828015453e+00;
									in->xdot[3] = 3.228986350219792e-04;
									in->xdot[4] = -2.067303680241406e+01;
									in->xdot[5] = -2.005255049361735e+02;
								}
								break;
							}
							case GNSF:  // GNSF
							{
								// set model funtions
								sim_set_model(config, in, "phi_fun", &phi_fun);
								sim_set_model(config, in, "phi_fun_jac_y", &phi_fun_jac_y);
								sim_set_model(config, in, "phi_jac_y_uhat", &phi_jac_y_uhat);
								sim_set_model(config, in, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);

								// import model matrices
								external_function_generic *get_model_matrices =
										(external_function_generic *) &get_matrices_fun;
								gnsf_model *model = (gnsf_model *) in->model;
								sim_gnsf_import_matrices(gnsf_dim, model, get_model_matrices);
								break;
							}
							// case NEW_LIFTED_IRK:  // new_lifted_irk
							// {
							//     sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
							//     sim_set_model(config, in, "impl_ode_fun_jac_x_xdot_u",
							//              &impl_ode_fun_jac_x_xdot_u);
							//     break;
							// }
							default :
							{
								printf("\nnot enough sim solvers implemented!\n");
								exit(1);
							}
						}

					/* seeds */
						for (int ii = 0; ii < nx * NF; ii++)
							in->S_forw[ii] = 0.0;
						for (int ii = 0; ii < nx; ii++)
							in->S_forw[ii * (nx + 1)] = 1.0;

						// seeds adj
						for (int ii = 0; ii < nx; ii++)
							in->S_adj[ii] = 1.0;
						for (int ii = nx; ii < nx + nu; ii++)
							in->S_adj[ii] = 0.0;
					/* sim solver  */
						sim_solver = sim_create(config, dims, opts);
						int acados_return;

						if (plan.sim_solver == GNSF){  // for gnsf: perform precomputation
							gnsf_model *model = (gnsf_model *) in->model;
							sim_gnsf_precompute(config, gnsf_dim, model, opts,
										sim_solver->mem, sim_solver->work, in->T);
						}

						
					/* sim loop */
						for (int ii=0; ii<nsim; ii++)
						{
							// update initial state
							for (int jj = 0; jj < nx; jj++)
								in->x[jj] = x_sim[ii*nx+jj];

							in->u[1] = 0.0;

							// update parameters
							switch (plan.sim_solver)
							{
								case ERK:  // ERK
								{
									expl_ode_fun.set_param(&expl_ode_fun, p_sim+ii*np);
									expl_vde_for.set_param(&expl_vde_for, p_sim+ii*np);
									expl_vde_for.set_param(&expl_vde_adj, p_sim+ii*np);
									break;
								}
								case IRK:  // IRK
								{
									impl_ode_fun.set_param(&impl_ode_fun, p_sim+ii*np);
									impl_ode_fun_jac_x_xdot.set_param(&impl_ode_fun_jac_x_xdot, p_sim+ii*np);
									impl_ode_jac_x_xdot_u.set_param(&impl_ode_jac_x_xdot_u, p_sim+ii*np);
									break;
								}
								case GNSF:  // GNSF
								{
									phi_fun.set_param(&phi_fun, p_sim+ii*np);
									phi_fun_jac_y.set_param(&phi_fun_jac_y, p_sim+ii*np);
									phi_jac_y_uhat.set_param(&phi_jac_y_uhat, p_sim+ii*np);
									f_lo_fun_jac_x1k1uz.set_param(&f_lo_fun_jac_x1k1uz, p_sim+ii*np);
									break;
								}
								default :
								{
									printf("\nnot enough sim solvers implemented!\n");
									exit(1);
								}
							}

							// execute simulation step with current input and state
							acados_return = sim_solve(sim_solver, in, out);
							if (acados_return != 0)
							{
								printf("error in sim solver\n");
								return ACADOS_FAILURE;
							}

							cpu_time[ii] = out->info->CPUtime;
							lss_time[ii] = out->info->LAtime;
							ad_time[ii]  = out->info->ADtime;
							la_time[ii]  = out->info->CPUtime - out->info->ADtime;;

							// extract state at next time step
							for (int jj = 0; jj < nx; jj++)
								x_sim[(ii+1)*nx+jj] = out->xn[jj];

							/* compute errors w.r.t. reference solution */
							if (ii == 0){
								// error sim
								for (int jj = 0; jj < nx; jj++){
									error[jj] = fabs(out->xn[jj] - x_ref_sol[jj]);
								}
								norm_error = onenorm(nx, 1, error);
								rel_error_x = norm_error / norm_x_ref;

								if ( opts->sens_forw ){     // error_S_forw
									norm_error_forw = 0.0;
									for (int jj = 0; jj < nx*NF; jj++){
										error_S_forw[jj] = fabs(S_forw_ref_sol[jj] - out->S_forw[jj]);
									}
									norm_error_forw = onenorm(nx, nx + nu, error_S_forw);
									rel_error_forw = norm_error_forw / norm_S_forw_ref;
								}


								if ( opts->sens_adj ){               // error_S_adj
									for (int jj = 0; jj < nx + nu; jj++){
										error_S_adj[jj] = S_adj_ref_sol[jj] - out->S_adj[jj];
									}
									norm_error_adj = onenorm(1, nx +nu, error_S_adj);
									rel_error_adj = norm_error_adj / norm_S_adj_ref;
								}

								if ( opts->output_z ){      // error_z
									for (int jj = 0; jj < nz; jj++){
										error_z[jj] = fabs(out->zn[jj] - z_ref_sol[jj]);
									}
									norm_error_z = onenorm(nz, 1, error_z);
									rel_error_z = norm_error_z / norm_z_ref;
								}

								if ( opts->sens_algebraic ){        // error_S_alg
									for (int jj = 0; jj < nz * (nx + nu); jj++){
										error_S_alg[jj] = fabs(out->S_algebraic[jj] - S_alg_ref_sol[jj]);
									}
									norm_error_sens_alg = onenorm(nz, nx + nu, error_S_alg);
									rel_error_alg = norm_error_sens_alg / norm_S_alg_ref;
								}
							}
						}
						// take average of nsim timings
						cpu_time_experiment = average_of_doubles(cpu_time, nsim);
						lss_time_experiment = average_of_doubles(lss_time, nsim);
						ad_time_experiment  = average_of_doubles(ad_time, nsim); 
						la_time_experiment  = average_of_doubles(la_time, nsim);

						/* printing - OFF */
						#if 0
							printf("\nxn: \n");
							d_print_exp_mat(1, nx, &x_sim[nsim*nx], 1);

							double *S_forw_out = NULL;
							if(opts->sens_forw){
								S_forw_out = out->S_forw;
								printf("\nS_forw_out: \n");
								d_print_exp_mat(nx, NF, S_forw_out, nx);
							}
							// printf("time split: %f ms CPU, %f ms LA, %f ms AD\n\n", cpu_time, lss_time, ad_time);
							printf("\naverage time for 1 simulation step: %f ms (AD time: %f ms (%5.2f%%))\n", 1e3*cpu_time_experiment, 1e3*ad_time_experiment, 1e2*ad_time_experiment/cpu_time_experiment);
							printf("time spent in integrator outside of casADi %f \n", 1e3*(cpu_time_experiment-ad_time_experiment));
						#endif



					/* store experiment results in array entries */
						// int i_experiment = (num_stages - min_num_stages) * steps_in_experiment + (num_steps - min_num_steps);
						int i_experiment = ((num_stages - min_num_stages) * steps_in_experiment + (tested_steps)) *
											newton_in_experiment + (newton_iter - min_newton); // maybe easier to increment :P
						// options
						experiment_solver[i_experiment] 		 	= nss;
						experiment_num_stages[i_experiment] 	 	= num_stages;
						experiment_num_steps[i_experiment] 	 		= num_steps;
						experiment_newton_it[i_experiment] 		 	= newton_iter;
						experiment_jac_reuse[i_experiment] 		 	= (double) jac_reuse;
						experiment_sens_forw[i_experiment] 		 	= (double) sens_forw;
						experiment_sens_adj[i_experiment] 		 	= (double) sens_adj;
						experiment_output_z[i_experiment] 			= (double) output_z;
						experiment_sens_alg[i_experiment] 		 	= (double) sens_alg;
						// errors
						experiment_error_sim[i_experiment]       	= rel_error_x;
						experiment_error_sens_forw[i_experiment] 	= rel_error_forw;
						experiment_error_sens_adj[i_experiment] 	= rel_error_adj;
						experiment_error_z[i_experiment] 			= rel_error_z;
						experiment_error_sens_alg[i_experiment] 	= rel_error_alg;
						// timings
						/* take minimal timing over n_executions */

						if ( cpu_time_experiment <=  experiment_cpu_time[i_experiment] ){
							experiment_cpu_time[i_experiment]  	= cpu_time_experiment;
							experiment_ad_time[i_experiment]   	= ad_time_experiment;
							experiment_lss_time[i_experiment]   = lss_time_experiment;
							experiment_la_time[i_experiment]  	= la_time_experiment;
						}
						// experiment_cpu_time[i_experiment]  	= experiment_cpu_time[i_experiment] < cpu_time_experiment ?
						// 				experiment_cpu_time[i_experiment] : cpu_time_experiment;
						// experiment_ad_time[i_experiment]   	= experiment_ad_time[i_experiment] < ad_time_experiment ?
						// 										experiment_ad_time[i_experiment] : ad_time_experiment;
						// experiment_lss_time[i_experiment]   = experiment_lss_time[i_experiment] < lss_time_experiment ?
						// 										experiment_lss_time[i_experiment] : lss_time_experiment;
						// experiment_la_time[i_experiment]  	= experiment_la_time[i_experiment] < la_time_experiment ?
						// 										experiment_la_time[i_experiment] : la_time_experiment;
					/* free memory */
						free(dims);
						free(sim_solver);
						free(in);
						free(out);
						free(opts);
						free(config);
				}  // end newton loop
			}  // end num_steps loop
		}  // end num_stages loop
	}  //  end n_execution loop

	/* print results to file */
		char export_filename[150] = "/home/oj/Git/1Thesis/1Matlab_prototypes/evaluation/results/results_";
		// append model name
		strcat(export_filename, "wt_nx6_");
		if (nss == 2){
			strcat(export_filename, "irk");
		}
		else if (nss == 3){
			strcat(export_filename, "gnsf");
		}
		else if (nss == 1){
			strcat(export_filename, "erk");
		}
		// append date identifier
		strcat(export_filename, "_september_3_1");
		// append additional identifier
		if (gnsf_init){
			strcat(export_filename, "_init_eq");
		}
		else {
			strcat(export_filename, "_init0");
		}
		// append file format
		strcat(export_filename, ".txt");

		// open file
		FILE *file_handle = fopen(export_filename, "wr");
		assert(file_handle != NULL);

		// 1st line -- options
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_solver	  , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_num_stages, 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_num_steps , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_newton_it , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_jac_reuse , 1);
		// 6th
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_sens_forw , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_sens_adj  , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_output_z  , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_sens_alg  , 1);
		/* errors */ 
		// line 10
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_error_sim ,  1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_error_sens_forw, 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_error_sens_adj, 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_error_z, 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_error_sens_alg, 1);
		/* timings */
		// line 15
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_cpu_time  , 1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_ad_time,    1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_lss_time   ,1);
		d_print_to_file_exp_mat(file_handle, 1, num_experiments, experiment_la_time ,   1);
		// line 19

		// close file
		fclose(file_handle);

	}  // end for solver loop

/* free functions */
	// explicit model
	external_function_param_casadi_free(&expl_ode_fun);
	external_function_param_casadi_free(&expl_vde_for);
	external_function_param_casadi_free(&expl_vde_adj);
	// implicit model
	external_function_param_casadi_free(&impl_ode_fun);
	external_function_param_casadi_free(&impl_ode_fun_jac_x_xdot);
	external_function_param_casadi_free(&impl_ode_jac_x_xdot_u);
	// gnsf model
	external_function_param_casadi_free(&f_lo_fun_jac_x1k1uz);
	external_function_param_casadi_free(&phi_fun);
	external_function_param_casadi_free(&phi_fun_jac_y);
	external_function_param_casadi_free(&phi_jac_y_uhat);
	external_function_casadi_free(&get_matrices_fun);

	double test_time = acados_toc(&test_timer);
	printf("\n********   INTEGRATOR TEST COMPLETE    ******* \ntotal runtime = %f  [min]\n\n", test_time/60);

	free(x_sim);

    return 0;
}
