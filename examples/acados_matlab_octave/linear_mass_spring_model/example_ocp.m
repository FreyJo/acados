%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

clear all;clc;
check_acados_requirements()

%% acados ocp model
model = get_linear_mass_spring_model();
nx = length(model.x);
nu = length(model.u);
ny = nx + nu;
ny_e = nx;

nbx = nx/2;
nbu = nu;
%% set up OCP
ocp = AcadosOcp();
ocp.model = model;

T = 10.0; % horizon length time
N = 20;

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'; % 'EXACT', 'GAUSS_NEWTON'
ocp.solver_options.regularize_method = 'NO_REGULARIZE';
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.sim_method_num_stages = 4;
ocp.solver_options.sim_method_num_steps = 3;

% cost
cost_type = 'AUTO';
ocp.cost.cost_type_0 = cost_type;
ocp.cost.cost_type = cost_type;
ocp.cost.cost_type_e = cost_type;

Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = eye(ny); for ii=1:nu W(ii,ii)=2.0; end % weight matrix in lagrange term
W_e = eye(ny_e); % weight matrix in mayer term
yr = zeros(ny, 1); % output reference in lagrange term
yr_e = zeros(ny_e, 1); % output reference in mayer term

sym_x = model.x;
sym_u = model.u;
yr_u = zeros(nu, 1);
yr_x = zeros(nx, 1);
dWu = 2*ones(nu, 1);
dWx = ones(nx, 1);
ymyr_0 = sym_u - yr_u;
ymyr = [sym_u; sym_x] - [yr_u; yr_x];
ymyr_e = sym_x - yr_x;

if (strcmp(cost_type, 'LINEAR_LS'))
    ocp.cost.Vu_0 = Vu;
    ocp.cost.Vu = Vu;
    ocp.cost.Vx_0 = Vx;
    ocp.cost.Vx = Vx;
    ocp.cost.Vx_e = Vx_e;
    ocp.cost.W_0 = W;
    ocp.cost.W = W;
    ocp.cost.W_e = W_e;
    ocp.cost.yref_0 = yr;
    ocp.cost.yref = yr; 
    ocp.cost.yref_e = yr_e;
elseif strcmp(cost_type, 'NONLINEAR_LS')
    ocp.model.cost_y_expr_0 = sym_u;
    ocp.model.cost_y_expr = [sym_u; sym_x];
    ocp.model.cost_y_expr_e = sym_x;
    ocp.cost.W_0 = W;
    ocp.cost.W = W;
    ocp.cost.W_e = W_e;
    ocp.cost.yref_0 = yr;
    ocp.cost.yref = yr;
    ocp.cost.yref_e = yr_e;
else
    cost_expr_ext_cost_0 = 0.5 * ymyr_0' * (dWu .* ymyr_0);
    cost_expr_ext_cost = 0.5 * ymyr' * ([dWu; dWx] .* ymyr);
    cost_expr_ext_cost_e = 0.5 * ymyr_e' * (dWx .* ymyr_e);

    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
    ocp.model.cost_expr_ext_cost = cost_expr_ext_cost;
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end

% constraints
x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;
ocp.constraints.x0 = x0;
ocp.constraints.idxbx = (0:nbx-1)';
ocp.constraints.lbx = -4*ones(nbx, 1);
ocp.constraints.ubx = 4*ones(nbx, 1);
ocp.constraints.idxbu = (0:nbu-1)';
ocp.constraints.lbu = -0.5*ones(nbu, 1);
ocp.constraints.ubu = 0.5*ones(nbu, 1);
%% acados ocp solver
% create ocp
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

%x0(1) = 1.5;
%ocp_solver.set('constr_x0', x0);

% if not set, the trajectory is initialized with the previous solution

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

if status==0
	fprintf('\nsuccess!\n\n');
else
	fprintf('\nsolution failed!\n\n');
end

% plot result
figure()
subplot(2, 1, 1)
plot(0:N, x);
title('trajectory')
ylabel('x')
subplot(2, 1, 2)
plot(1:N, u);
ylabel('u')
xlabel('sample')

if is_octave()
    waitforbuttonpress;
end
