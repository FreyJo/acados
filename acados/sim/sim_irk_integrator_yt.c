// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_irk_integrator_yt.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"


int irk_calculate_memory_size(sim_RK_opts *opts, sim_in *in)
{

    int nx = in->nx;
    int nu = in->nu; 
    // int NF = in->NF;
    
    int ns = opts->num_stages; // number of stages

    int size = sizeof(sim_irk_memory);

    size += 5*sizeof(struct d_strmat); // JG, JGf, JKf, JFK
    size += 4*sizeof(struct d_strvec); // rG, K, xt, xn
    size += 2*sizeof(struct d_strvec); // lambda,lambdaK

    size += d_size_strmat(nx*ns, nx*ns); // JG
    size += 2*d_size_strmat(nx*ns, nx+nu); // JGf, JKf
    size += d_size_strmat(nx, nx*ns); // JFK
    size += d_size_strmat(nx, nx+nu); // S_forw

    size += 4*d_size_strvec(nx*ns); // rG, K
    size += 2*d_size_strvec(nx); // xt, x

    size += d_size_strvec(nx); // lambda
    size += d_size_strvec(nx*ns); // lambdaK

    size += nx * sizeof(double); //  rGt
    size += nx * (2*nx+nu) * sizeof(double); // jac_out
    size += nx * nx * sizeof(double); // Jt
    size += (2*nx + nu) * sizeof(double); // ode_args

    size += nx *ns * sizeof(int); // ipiv

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}

char *assign_irk_memory(sim_RK_opts *opts, sim_in *in, sim_irk_memory **memory, void *raw_memory)
{

    int nx = in->nx;
    int nu = in->nu; 
    // int NF = in->NF;
        
    int ns = opts->num_stages; // number of stages

    char *c_ptr = (char *)raw_memory;

    *memory = (sim_irk_memory *) c_ptr;
    c_ptr += sizeof(sim_irk_memory);

    (*memory)->JG = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);   

    (*memory)->JGf = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->JKf = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->JFK = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->S_forw = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->rG = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    (*memory)->K = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    (*memory)->xt = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

     (*memory)->xn = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

     (*memory)->lambda = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    (*memory)->lambdaK = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    struct d_strmat *JG = (*memory)->JG; 
    struct d_strmat *JGf = (*memory)->JGf;
    struct d_strmat *JKf = (*memory)->JKf;
    struct d_strmat *JFK = (*memory)->JFK;
    struct d_strmat *S_forw = (*memory)->S_forw;

    struct d_strvec *rG = (*memory)->rG;
    struct d_strvec *K = (*memory)->K;
    struct d_strvec *xt = (*memory)->xt;
    struct d_strvec *xn = (*memory)->xn;

    struct d_strvec *lambda = (*memory)->lambda;
    struct d_strvec *lambdaK = (*memory)->lambdaK;
    

    // align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    d_create_strmat(nx*ns, nx*ns, JG, c_ptr);
    c_ptr += d_size_strmat(nx*ns, nx*ns);

    d_create_strmat(nx*ns, nx+nu, JGf, c_ptr);
    c_ptr += d_size_strmat(nx*ns, nx+nu);

    d_create_strmat(nx*ns, nx+nu, JKf, c_ptr);
    c_ptr += d_size_strmat(nx*ns, nx+nu);

    d_create_strmat(nx, nx*ns, JFK, c_ptr);
    c_ptr += d_size_strmat(nx, nx*ns);

    d_create_strmat(nx, nx+nu, S_forw, c_ptr);
    c_ptr += d_size_strmat(nx, nx+nu);

    d_create_strvec(nx*ns, rG, c_ptr);
    c_ptr += d_size_strvec(nx*ns);

    d_create_strvec(nx*ns, K, c_ptr);
    c_ptr += d_size_strvec(nx*ns);

    d_create_strvec(nx, xt, c_ptr);
    c_ptr += d_size_strvec(nx);

    d_create_strvec(nx, xn, c_ptr);
    c_ptr += d_size_strvec(nx);

    d_create_strvec(nx, lambda, c_ptr);
    c_ptr += d_size_strvec(nx);

    d_create_strvec(nx, lambdaK, c_ptr);
    c_ptr += d_size_strvec(nx*ns);

    (*memory)->rGt = (double *)c_ptr;
    c_ptr += nx * sizeof(double);

    (*memory)->jac_out = (double *)c_ptr;
    c_ptr += nx * (2*nx+nu) * sizeof(double);

    (*memory)->Jt = (double *)c_ptr;
    c_ptr += nx * nx * sizeof(double);

    (*memory)->ode_args = (double *)c_ptr;
    c_ptr += (2*nx + nu) * sizeof(double);

    (*memory)->ipiv = (int *)c_ptr;
    c_ptr += nx * ns * sizeof(int);

    return c_ptr;
}

sim_irk_memory *sim_irk_create_memory(sim_RK_opts *opts, sim_in *in)
{
    sim_irk_memory *memory;

    int_t bytes = irk_calculate_memory_size(opts, in);
    void *ptr = malloc(bytes);
    char *ptr_end = assign_irk_memory(opts, in, &memory, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    return memory;
}


int sim_irk_yt(const sim_in *in, sim_out *out, void *opts_, void *mem_){
    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    int ii, jj, iter, kk, ss;
    double a,b;

    int nx = in->nx;
	int nu = in->nu;
	// int NF = in->NF;
    int num_steps = in->num_steps;
    double *x = in->x;
	double *u = in->u;
    double step = in->step;
    double *S_forw_in = in->S_forw;
    double *S_adj_in = in->S_adj;
    
    int num_stages = opts->num_stages;
    int newton_iter = opts->newton_iter;
    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
 
	double *rGt = mem->rGt;  // residual of one stage of the impl. ODE
	double *jac_out = mem->jac_out;
    double *Jt = mem->Jt;
    double *ode_args = mem->ode_args;
    int *ipiv = mem->ipiv;
	struct d_strmat *JG = mem->JG;
	struct d_strvec *rG = mem->rG; // residual of the of the IRK equations G
    struct d_strvec *K = mem->K;
    struct d_strmat *JGf = mem->JGf;
    struct d_strmat *JKf = mem->JKf;
    struct d_strvec *xt = mem->xt;
    struct d_strvec *xn = mem->xn;
    struct d_strmat *JFK = mem->JFK;
    struct d_strmat *S_forw = mem->S_forw;

    struct d_strvec *lambda = mem->lambda;
    struct d_strvec *lambdaK = mem->lambdaK;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;

    acados_timer timer;

    // initialize
    dgese_libstr(nx*num_stages, nx*num_stages, 0.0, JG, 0, 0);
    dgese_libstr(nx*num_stages, nx+nu, 0.0, JGf, 0, 0);
    dgese_libstr(nx*num_stages, nx+nu, 0.0, JKf, 0, 0);
    for (ii=0;ii<num_stages;ii++){
        b = step * b_vec[ii];
        ddiare_libstr(nx, b, JFK, 0, ii*nx);
    }
    d_cvt_mat2strmat(nx, nx+nu, S_forw_in, nx, S_forw, 0, 0);

    dvecse_libstr(nx*num_stages, 0.0, rG, 0);
    dvecse_libstr(nx*num_stages, 0.0, K, 0);
    d_cvt_vec2strvec(nx, x, xn, 0);

    dvecse_libstr(nx*num_stages, 0.0, lambdaK, 0);
    d_cvt_vec2strvec(nx, S_adj_in, lambda, 0);
    
    for (kk=0;kk<2*nx;kk++) //initialize x,xdot with zeros
        ode_args[kk] = 0.0;
    for (kk=0;kk<nu;kk++) //set controls
        ode_args[2*nx+kk] = u[kk];

    for (kk=0;kk<nx;kk++){
        rGt[kk] = 0.0;
    }
    for (kk=0;kk<nx*(2*nx+nu);kk++)
        jac_out[kk] = 0.0;
    for (kk=0;kk<nx*nx;kk++)
        Jt[kk] = 0.0;
    for (kk=0;kk<num_stages*nx;kk++)
        ipiv[kk] = kk; //ipiv?
        
    // start the loop
    acados_tic(&timer);
    for(ss=0; ss<num_steps; ss++){

        //  obtain Kn
        for(iter=0; iter<newton_iter; iter++){
           
            for(ii=0; ii<num_stages; ii++){ // ii-th row of tableau   
                
                // take x(n); copy a strvec into a strvec
                dveccp_libstr(nx, xn, 0, xt, 0);

                for(jj=0; jj<num_stages; jj++){ // jj-th col of tableau
                    a = A_mat[ii+num_stages*jj];
					if(a!=0){           // xt = xt + T_int * a[i,j]*K_j
                        a *= step; 
                        daxpy_libstr(nx, a, K, jj*nx, xt, 0, xt, 0); 
					}
                }
                // put xn+sum kj into first nx elements of ode_arg               
                d_cvt_strvec2vec(nx, xt, 0, ode_args);
                
                // put ki into the next nx elements of ode_args
                d_cvt_strvec2vec(nx, K, ii*nx, ode_args+nx);      

                // compute the residual of implicit ode at time t_ii, store value in rGt 
                in->eval_impl_res(nx, nu, ode_args, rGt, in->impl_ode);
                // fill in elements of rG  - store values rGt on (ii*nx)th position of rG
                d_cvt_vec2strvec(nx, rGt, rG, ii*nx);
                
                // compute the jacobian of implicit ode
                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
                // compute the blocks of JG
                for (jj=0; jj<num_stages; jj++){ //compute the block (ii,jj)th block = Jt
                    a = A_mat[ii + num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0; kk<nx*nx; kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    // fill in the ii-th, jj-th block of JG
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);
                } // end jj

            } // end ii
            
            //DGETRF computes an LU factorization of a general M-by-N matrix A
            //using partial pivoting with row interchanges. // and store it in JG
            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv);
            // permute also the r.h.s
            dvecpe_libstr(nx*num_stages, ipiv, rG, 0);

            // solve JG * y = rG, JG on the (l)eft, (l)ower-trian, (n)o-trans
            //                    (u)nit trian
            dtrsv_lnu_libstr(nx*num_stages, JG, 0, 0, rG, 0, rG, 0);
            
            // solve JG * x = rG, JG on the (l)eft, (u)pper-trian, (n)o-trans
            //                    (n)o unit trian , and store x in rG
            dtrsv_unn_libstr(nx*num_stages, JG, 0, 0, rG, 0, rG, 0);
            // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is DeltaK
            daxpy_libstr(nx*num_stages, -1.0, rG, 0, K, 0, K, 0);
        }// end iter
        
        // evaluate forward sensitivities
        if (in->sens_forw){

            // evaluate JG(xn,Kn)         
            for(ii=0; ii<num_stages; ii++){  

                dveccp_libstr(nx, xn, 0, xt, 0);
                //compute xt = final x;
                for(jj=0; jj<num_stages; jj++){ 
                    a = A_mat[ii+num_stages*jj];
                    if(a!=0){
                        a *= step;
                        daxpy_libstr(nx, a, K, jj*nx, xt, 0, xt, 0);
                    }
                }               
                d_cvt_strvec2vec(nx, xt, 0, ode_args);                                           
                d_cvt_strvec2vec(nx, K, ii*nx, ode_args+nx);  

                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                d_cvt_mat2strmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);

                in->eval_impl_jac_u(nx, nu, ode_args, jac_out+2*nx*nx, in->impl_jac_u);
                d_cvt_mat2strmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);

                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
                
                // maybe this part can also be written using blasfeo
                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ii+num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);            
                } // end jj     
            } // end ii

            // factorize JG
            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv);

            // obtain JKf
            dgemm_nn_libstr(nx*num_stages, nx, nx, 1.0, JGf, 0, 0, S_forw, 0, 0, 0.0, JKf, 0, 0, JKf, 0, 0);
            dgemm_nn_libstr(nx*num_stages, nu, nx, 1.0, JGf, 0, 0, S_forw, 0, nx, 1.0, JGf, 0, nx, JKf, 0, nx);  
            drowpe_libstr(nx*num_stages, ipiv, JKf);
            dtrsm_llnu_libstr(nx*num_stages, nx+nu, 1.0, JG, 0, 0, JKf, 0, 0, JKf, 0, 0); 
            dtrsm_lunn_libstr(nx*num_stages, nx+nu, 1.0, JG, 0, 0, JKf, 0, 0, JKf, 0, 0);

            // update forward sensitivity
            dgemm_nn_libstr(nx, nx+nu, nx*num_stages, -1.0, JFK, 0, 0, JKf, 0, 0, 1.0, S_forw, 0, 0, S_forw, 0, 0);
        }

        if (in->sens_adj){        
            for(ii=0; ii<num_stages; ii++){  

                dveccp_libstr(nx, xn, 0, xt, 0);

                for(jj=0; jj<num_stages; jj++){ 
                    a = A_mat[ii+num_stages*jj];
                    if(a!=0){
                        a *= step;
                        daxpy_libstr(nx, a, K, jj*nx, xt, 0, xt, 0);
                    }
                }               
                d_cvt_strvec2vec(nx, xt, 0, ode_args);                                           
                d_cvt_strvec2vec(nx, K, ii*nx, ode_args+nx);  

                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                d_cvt_mat2strmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);

                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
            
                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ii+num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);            
                } // end jj     
            } // end ii

            // factorize JG
            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv);

            dgemv_t_libstr(nx, nx*num_stages, 1.0, JFK, 0, 0, lambda, 0, 0.0, lambda, 0, lambda, 0);

            dtrsv_lnn_libstr(nx*num_stages, JG, 0, 0, lambda, 0, lambda, 0);

            // neither is implemented
            dtrsv_unu_libstr(nx*num_stages, JG, 0, 0, lambda, 0, lambda, 0);
            dtrsv_ltu_libstr(nx*num_stages, JG, 0, 0, lambda, 0, lambda, 0);

            // dvecpei_libstr();


        }

        // obtain x(n+1)
        for(ii=0;ii<num_stages;ii++){
            b = step * b_vec[ii];
            daxpy_libstr(nx, b, K, ii*nx, xn, 0, xn, 0);          
        }
            
            
    }// end int step ss

    // extract output

    out->info->CPUtime = acados_toc(&timer);

    d_cvt_strvec2vec(nx, xn, 0, x_out);

    d_cvt_strmat2mat(nx, nx+nu, S_forw, 0, 0, S_forw_out, nx);

    return 0;
}