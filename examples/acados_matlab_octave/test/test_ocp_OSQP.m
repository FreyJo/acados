%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface

import casadi.*
addpath('../pendulum_on_cart_model/');


%% discretization
N = 20;
T = 1; % time horizon length
x0 = [0; pi; 0; 0];

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_osqp';
    % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases, partial_condensing_osqp
qp_solver_cond_N = 5; % for partial condensing
% integrator type
sim_method = 'erk'; % erk, irk, irk_gnsf

%% model dynamics
old_model = pendulum_on_cart_model();
nx = old_model.nx;
nu = old_model.nu;
model = AcadosModel();
model.name = 'pendulum';
model.x = old_model.sym_x;
model.xdot = old_model.sym_xdot;
model.u = old_model.sym_u;
model.cost_expr_ext_cost = old_model.cost_expr_ext_cost;
model.cost_expr_ext_cost_e = old_model.cost_expr_ext_cost_e;
model.con_h_expr = old_model.constr_expr_h;
model.con_h_expr_0 = old_model.constr_expr_h_0;
if strcmp(sim_method, 'erk')
    model.f_expl_expr = old_model.dyn_expr_f_expl;
else
    model.f_impl_expr = old_model.dyn_expr_f_impl;
end

%% acados OCP
ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = upper(nlp_solver);
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.integrator_type = upper(sim_method);
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.qp_solver_iter_max = 2000;
ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
ocp.code_gen_options.ext_fun_compile_flags = '';

U_max = 80;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.x0 = x0;

%% create ocp solver
ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);

%% call ocp solver
% update initial state
ocp_solver.set('constr_x0', x0);

% set trajectory initialization
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);
ocp_solver.set('init_pi', zeros(nx, N))

% change values for specific shooting node using:
%   ocp_solver.set('field', value, optional: stage_index)
ocp_solver.set('constr_lbx', x0, 0)

% solve
ocp_solver.solve();
% get solutionn
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

status = ocp_solver.get('status'); % 0 - success
ocp_solver.print('stat')

if status == 0
    disp('test_ocp_OSQP: success!');
else
    error(['test_ocp_OSQP: Failed with status ', num2str(status)]);
end
