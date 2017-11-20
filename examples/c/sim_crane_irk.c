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

#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"

int main() {
    int ii;
    int jj;
    
    int nx = 8;
    int nu = 2;
    int NF = nx + nu; // columns of forward seed

    real_t T = 0.1;
    real_t *xref;
    xref = (double*)calloc(nx, sizeof(real_t));
    xref[2] = 0.8;

    double u_0[] = {40.108149413030752, -50.446662212534974};

    sim_in *in = create_sim_in(nx, nu ,NF);
// why always twice; with(out) eval ?!
    in->num_steps = 1;
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
    in->sens_adj = false;

    for (ii = 0; ii < nx; ii++) {
        in->x[ii] = xref[ii];
    }
    
    for (ii = 0; ii < nu; ii++) {
        in->u[ii] = u_0[ii];
    }


    for (ii = 0; ii < nx * NF; ii++) //initialize with unit matrix
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    /*double A_impl[] = {0.1389, 0.3003, 0.2680 ,
        -0.0360, 0.2222, 0.4804,
        0.0098, -0.0225, 0.1389};
    double B_impl[] = {0.2778, 0.4444, 0.2778}; */
    // double C_impl[] = {0.0, 0.5, 0.5, 0.0};
    double A_impl[] = {   0.086961218376613,  -0.026602569005020,   0.012626331546231,  -0.003554980917825,
        0.188113390907761,   0.163029269279675,  -0.027877891817945,   0.006735231630509,
        0.167184974461456,   0.353952685725980,   0.163045524628359,  -0.014193184815796,
        0.177475187009790,   0.313448462361803,   0.352677362913056,   0.086958987715352};

    double B_impl[] = {0.173919066836852, 0.326078801640426,  0.326066657180823, 0.173935474341899};
    double C_impl[] = {   0.069430000000000,   0.330000000000000,   0.669990000000000,   0.930560000000000};

    int num_stages = 4;
    sim_RK_opts *irk_opts = create_sim_RK_opts(num_stages);
    irk_opts->newton_iter = 30;

    for (ii=0;ii<num_stages*num_stages;ii++){
        irk_opts->A_mat[ii] = A_impl[ii];
    }
    for (ii=0;ii<num_stages;ii++){
        irk_opts->b_vec[ii] = B_impl[ii];
        irk_opts->c_vec[ii] = C_impl[ii];
    }

    sim_irk_memory *irk_mem = sim_irk_create_memory(irk_opts, in);

    sim_out *out = create_sim_out(nx, nu, NF);

    int flag = sim_irk_yt(in, out, irk_opts, irk_mem);

    double *xn = out->xn;
    double *S_forw_out = out->S_forw;

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
    
    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime*1000);
    //printf("LAt: %8.4f [ms]\n", out->info->LAtime*1000);

    free(xref);
    free(in);
    free(irk_opts);
    free(out);
    free(irk_mem);
    
    return flag;
}