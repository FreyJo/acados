
%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface

addpath('../linear_mass_spring_model/');

%% arguments
N = 20;
shooting_nodes = [ linspace(0,1,N/2) linspace(1.1,5,N/2+1) ];

model_name = 'lin_mass';

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
% regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'mirror';
regularize_method = 'convexify';
nlp_solver_max_iter = 100;
nlp_solver_ext_qp_res = 1;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_cond_N = 5;
% sim_method = 'irk';
%sim_method = 'irk_gnsf';
sim_method = 'discrete';
sim_method_num_stages = 4 * ones(N,1);
sim_method_num_stages(end) = 2;
sim_method_num_steps = 3;
%cost_type = 'linear_ls';
%cost_type = 'nonlinear_ls';

%% create model entries
model = get_linear_mass_spring_model();
model.name = model_name;


% dims
T = 10.0; % horizon length time
nx = length(model.x);
nu = length(model.u);
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
nbx = 0;
nbu = 0;
ng = 0;
nh = nu+nx;
nh_e = nx;

% cost
Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = eye(ny); for ii=1:nu W(ii,ii)=2.0; end % weight matrix in lagrange term
W_e = eye(ny_e); % weight matrix in mayer term
yr = zeros(ny, 1); % output reference in lagrange term
yr_e = zeros(ny_e, 1); % output reference in mayer term
% constraints
x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;
if (ng>0)
	D = zeros(ng, nu); for ii=1:nu D(ii,ii)=1.0; end
	C = zeros(ng, nx); for ii=1:ng-nu C(nu+ii,ii)=1.0; end
	lg = zeros(ng, 1); for ii=1:nu lg(ii)=-0.5; end; for ii=1:ng-nu lg(nu+ii)=-4.0; end
	ug = zeros(ng, 1); for ii=1:nu ug(ii)= 0.5; end; for ii=1:ng-nu ug(nu+ii)= 4.0; end
	C_e = zeros(ng_e, nx); for ii=1:ng_e C_e(ii,ii)=1.0; end
	lg_e = zeros(ng_e, 1); for ii=1:ng_e lg_e(ii)=-4.0; end
	ug_e = zeros(ng_e, 1); for ii=1:ng_e ug_e(ii)= 4.0; end
elseif (nh>0)
	lh = zeros(nh, 1); for ii=1:nu lh(ii)=-0.5; end; for ii=1:nx lh(nu+ii)=-4.0; end
	uh = zeros(nh, 1); for ii=1:nu uh(ii)= 0.5; end; for ii=1:nx uh(nu+ii)= 4.0; end
	lh_e = zeros(nh_e, 1); for ii=1:nh_e lh_e(ii)=-4.0; end
	uh_e = zeros(nh_e, 1); for ii=1:nh_e uh_e(ii)= 4.0; end
else
	Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
	lbx = -4*ones(nbx, 1);
	ubx =  4*ones(nbx, 1);
	Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
	lbu = -0.5*ones(nbu, 1);
	ubu =  0.5*ones(nbu, 1);
end


%% acados OCP
ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.shooting_nodes = shooting_nodes;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(nlp_solver);
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.regularize_method = upper(regularize_method);
ocp.solver_options.nlp_solver_ext_qp_res = nlp_solver_ext_qp_res;
ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.integrator_type = upper(sim_method);
ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = sim_method_num_steps;
if strcmp(qp_solver, 'partial_condensing_hpipm')
	ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
end
ocp.code_gen_options.ext_fun_compile_flags = '';

ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.cost.cost_type_0 = 'EXTERNAL';
ocp.model.cost_expr_ext_cost_0 = 0.5 * model.u' * (2 * eye(nu)) * model.u;
ocp.model.cost_expr_ext_cost = 0.5 * ([model.u; model.x] .* [2*ones(nu,1); ones(nx,1)])' * ...
    ([model.u; model.x] .* [2*ones(nu,1); ones(nx,1)]);
ocp.model.cost_expr_ext_cost_e = 0.5 * model.x' * model.x;

ocp.model.con_h_expr = [model.u; model.x];
ocp.model.con_h_expr_e = model.x;
ocp.constraints.lh = lh;
ocp.constraints.uh = uh;
ocp.constraints.lh_e = lh_e;
ocp.constraints.uh_e = uh_e;
ocp.constraints.x0 = x0;

%% acados OCP solver
ocp_solver = AcadosOcpSolver(ocp);


% set trajectory initialization
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);


% solve
tic;
ocp_solver.solve();
time_ext = toc;

% store and load iterate
filename = 'iterate.json';
ocp_solver.store_iterate(filename, true);
ocp_solver.load_iterate(filename);
delete(filename)

% test QP dump
filename = 'qp.json';
ocp_solver.dump_last_qp_to_json(filename, false, 'Matlab')
delete(filename)

filename = 'qp_c_backend.json';
ocp_solver.dump_last_qp_to_json(filename, false)
delete(filename)

% test qp_diagnostics
qp_diagnostics_result = ocp_solver.qp_diagnostics();


% get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

% get info
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

% print statistics
ocp_solver.print('stat')

if status~=0
    error('ocp_nlp solver returned status nonzero');
elseif sqp_iter > 2
    error('ocp can be solved in 2 iterations!');
else
	fprintf(['\ntest_ocp_linear_mass_spring: success with sim method ', ...
        sim_method, ' !\n']);
end

% plot result
% figure()
% subplot(2, 1, 1)
% plot(0:N, x);
% title('trajectory')
% ylabel('x')
% subplot(2, 1, 2)
% plot(1:N, u);
% ylabel('u')
% xlabel('sample')
% 
% if is_octave()
%     waitforbuttonpress;
% end
