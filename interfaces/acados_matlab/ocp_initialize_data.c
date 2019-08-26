/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"
#include "mex_macros.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	long long *ptr;
	int ii, jj, kk, acados_size, matlab_size;
	mxArray *mex_field;
	char fun_name[50] = "ocp_initialize_data";
	char buffer [100]; // for error messages
    char field[30];

    const mxArray *matlab_model = prhs[0];
    const mxArray *matlab_options = prhs[1];
    const mxArray *C_struct = prhs[2];
    const mxArray *value_array;

	/* RHS */
	// config
	ptr = (long long *) mxGetData( mxGetField( C_struct, 0, "config" ) );
	ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
	// dims
	ptr = (long long *) mxGetData( mxGetField( C_struct, 0, "dims" ) );
	ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
	// opts
	ptr = (long long *) mxGetData( mxGetField( C_struct, 0, "opts" ) );
	void *opts = (void *) ptr[0];
	// in
	ptr = (long long *) mxGetData( mxGetField( C_struct, 0, "in" ) );
	ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
	// out
	ptr = (long long *) mxGetData( mxGetField( C_struct, 0, "out" ) );
	ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];

    // dims
    int ns = ocp_nlp_dims_get_from_attr(config, dims, out, 1, "ns");
	int N = dims->N;


    /* slacks */
    d_ptr = malloc(ns*sizeof(double));
	// Z
    field = "cost_Z";
    value_array = mxGetField( matlab_model, 0, field );
	if(value_array!=NULL)
    {
        matlab_size = (int) mxGetNumberOfElements( value_array );
        acados_size = ns * N;
        MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);

		value = mxGetPr( value_array );
		for(ii=0; ii<ns; ii++)
        {
			d_ptr[ii] = value[ii+ns*ii];
        }
		for(ii=0; ii<N; ii++)
        {
			ocp_nlp_cost_model_set(config, dims, in, ii, "Z", d_ptr);
        }
    }


	// Z_e
    field = "cost_Z_e";
    value_array = mxGetField( matlab_model, 0, field );
	if(value_array!=NULL)
    {
        matlab_size = (int) mxGetNumberOfElements( value_array );
        acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, N, "ns");
        MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);

		value = mxGetPr( value_array );

		for(ii=0; ii<ns; ii++)
        {
			d_ptr[ii] = value[ii+ns*ii];
        }
		ocp_nlp_cost_model_set(config, dims, in, N, "Z", d_ptr);
    }




	// Zl
    field = "cost_Zl";
    value_array = mxGetField( matlab_model, 0, field );
	if(value_array!=NULL)
    {
        matlab_size = (int) mxGetNumberOfElements( value_array );
        acados_size = N*;
        MEX_DIM_CHECK(fun_name, field, matlab_size, acados_size);

		value = mxGetPr( value_array );

		for(ii=0; ii<ns; ii++)
        {
			d_ptr[ii] = Zl[ii+ns*ii];
        }
		for(ii=0; ii<N; ii++)
        {
			ocp_nlp_cost_model_set(config, dims, in, ii, "Zl", d_ptr);
        }
		ocp_nlp_cost_model_set(config, dims, in, N, "Z", d_ptr);
    }

	if(set_Zl)
		{
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zl[ii+ns*ii];
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Zl", d_ptr);
			}
		free(d_ptr);
		}

	// Zl
	if(mxGetField( matlab_model, 0, "cost_Zl" )!=NULL)
		{
		set_Zl = true;
		Zl = mxGetPr( mxGetField( matlab_model, 0, "cost_Zl" ) );
		}
	// Zl_e
	if(mxGetField( matlab_model, 0, "cost_Zl_e" )!=NULL)
		{
		set_Zl_e = true;
		Zl_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Zl_e" ) );
		}
	// Zu
	if(mxGetField( matlab_model, 0, "cost_Zu" )!=NULL)
		{
		set_Zu = true;
		Zu = mxGetPr( mxGetField( matlab_model, 0, "cost_Zu" ) );
		}
	// Zu_e
	if(mxGetField( matlab_model, 0, "cost_Zu_e" )!=NULL)
		{
		set_Zu_e = true;
		Zu_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Zu_e" ) );
		}
	// z
	if(mxGetField( matlab_model, 0, "cost_z" )!=NULL)
		{
		set_z = true;
		z = mxGetPr( mxGetField( matlab_model, 0, "cost_z" ) );
		}
	// z_e
	if(mxGetField( matlab_model, 0, "cost_z_e" )!=NULL)
		{
		set_z_e = true;
		z_e = mxGetPr( mxGetField( matlab_model, 0, "cost_z_e" ) );
		}
	// zl
	if(mxGetField( matlab_model, 0, "cost_zl" )!=NULL)
		{
		set_zl = true;
		zl = mxGetPr( mxGetField( matlab_model, 0, "cost_zl" ) );
		}
	// zl_e
	if(mxGetField( matlab_model, 0, "cost_zl_e" )!=NULL)
		{
		set_zl_e = true;
		zl_e = mxGetPr( mxGetField( matlab_model, 0, "cost_zl_e" ) );
		}
	// zu
	if(mxGetField( matlab_model, 0, "cost_zu" )!=NULL)
		{
		set_zu = true;
		zu = mxGetPr( mxGetField( matlab_model, 0, "cost_zu" ) );
		}
	// zu_e
	if(mxGetField( matlab_model, 0, "cost_zu_e" )!=NULL)
		{
		set_zu_e = true;
		zu_e = mxGetPr( mxGetField( matlab_model, 0, "cost_zu_e" ) );
		}






	if(set_Zl_e)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zl_e[ii+ns*ii];
			}
		ocp_nlp_cost_model_set(config, dims, in, N, "Zl", d_ptr);
		free(d_ptr);
		}
	if(set_Zu)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zu[ii+ns*ii];
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Zu", d_ptr);
			}
		free(d_ptr);
		}
	if(set_Zu_e)
		{
		d_ptr = malloc(ns*sizeof(double));
		for(ii=0; ii<ns; ii++)
			{
			d_ptr[ii] = Zu_e[ii+ns*ii];
			}
		ocp_nlp_cost_model_set(config, dims, in, N, "Zu", d_ptr);
		free(d_ptr);
		}
	if(set_z)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "z", z);
			}
		}
	if(set_z_e)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "z", z_e);
		}
	if(set_zl)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "zl", zl);
			}
		}
	if(set_zl_e)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "zl", zl_e);
		}
	if(set_zu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "zu", zu);
			}
		}
	if(set_zu_e)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "zu", zu_e);
		}


    free(d_ptr);
}