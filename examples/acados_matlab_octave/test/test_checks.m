%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

%% test of native matlab interface
clear VARIABLES

addpath('../linear_mass_spring_model/');

%% arguments
method = 'irk';
sens_forw = true;
num_stages = 4;
num_steps = 4;

Ts = 0.1;

%% model
model = get_linear_mass_spring_model();

model_name = ['lin_mass_' method];
nx = length(model.x);
nu = length(model.u);

%% acados sim
sim = AcadosSim();
sim.model = model;
sim.model.name = model_name;
sim.solver_options.Tsim = Ts;
sim.solver_options.integrator_type = upper(method);
sim.solver_options.num_stages = num_stages;
sim.solver_options.num_steps = num_steps;
sim.solver_options.sens_forw = sens_forw;
sim.model.f_impl_expr = model.f_impl_expr;
sim_solver = AcadosSimSolver(sim);

% Note: this does not work with gnsf, because it needs to be available
% in the precomputation phase
% 	sim_solver.set('T', Ts);

%% test check, this should fail!
try
    sim_solver.set('x', zeros(nx+1, 1));
    error('test_checks: setter accepted a state with the wrong dimension');
catch exception
    if contains(exception.message, 'setter accepted')
        rethrow(exception);
    end
end
