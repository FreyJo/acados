#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_erk_integrator_yt.h"
#include "acados/sim/sim_irk_integrator_yt.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/yutao_model/yutao_model.h"

// blasfeo
#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"


#define M_PI 3.14159265358979323846

int main() {
    int ii;
    int jj;
    
    int nx = 4;
    int nu = 1;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

    sim_in *in = create_sim_in(nx, nu ,NF);

    in->num_steps = 4;
    in->step = T / in->num_steps;
    in->impl_ode = &impl_odeFun;
    in->eval_impl_res = &impl_ode_fun;
    in->impl_jac_x = &impl_jacFun_x;
    in->eval_impl_jac_x = &impl_jac_x_fun;
    in->impl_jac_xdot = &impl_jacFun_xdot;
    in->eval_impl_jac_xdot = &impl_jac_xdot_fun;
    in->impl_jac_u = &impl_jacFun_u;
    in->eval_impl_jac_u = &impl_jac_u_fun;

    in->sens_forw = true;
    in->sens_adj = true;

    for (ii = 0; ii < nx; ii++) {
        in->x[ii] = xref[ii];
    }
    for (ii = 0;ii < nu; ii++){
        in->u[ii] = 1.0;
    }

    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    for (ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;
    // for (ii = 0; ii < nu; ii++)
    //     in->S_adj[nx+ii] = 0.0;

    double A_impl[] = {0.1389, 0.3003, 0.2680 ,
        -0.0360, 0.2222, 0.4804,
        0.0098, -0.0225, 0.1389};
    double B_impl[] = {0.2778, 0.4444, 0.2778};
    // double C_impl[] = {0.0, 0.5, 0.5, 0.0};

    int num_stages = 3;
    sim_RK_opts *irk_opts = create_sim_RK_opts(num_stages);
    irk_opts->newton_iter = 3;

    for (ii=0;ii<num_stages*num_stages;ii++){
        irk_opts->A_mat[ii] = A_impl[ii];
    }
    for (ii=0;ii<num_stages;ii++){
        irk_opts->b_vec[ii] = B_impl[ii];
        // irk_opts->c_vec[ii] = C_impl[ii];
    }

    sim_irk_memory *irk_mem = sim_irk_create_memory(irk_opts, in);

    sim_out *out = create_sim_out(nx, nu, NF);

    int flag = sim_irk_yt(in, out, irk_opts, irk_mem);

    double *xn = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    printf("\nS_forw_out: \n");
    for (ii=0;ii<nx;ii++){
        for (jj=0;jj<NF;jj++)
            printf("%8.5f ",S_forw_out[jj*nx+ii]);
        printf("\n");
    }

    printf("\nS_adj_out: \n");
    for (ii=0;ii<nx+nu;ii++){
        printf("%8.5f ",S_adj_out[ii]);
    }
    printf("\n"); 
    
    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime*1000);

    free(xref);
    free(in);
    free(irk_opts);
    free(out);
    free(irk_mem);
    
    return flag;
}