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

for integrator = {'irk_gnsf', 'irk', 'erk'}
    method = integrator{1};

    %% arguments
    sens_forw = true;
    jac_reuse = true;
    num_stages = 3;
    num_steps = 4;
    newton_iter = 3;
    gnsf_detect_struct = 'true';

    Ts = 0.1;
    FD_epsilon = 1e-6;

    %% model
    model = linear_mass_spring_model();

    model_name = ['lin_mass_' method];
    nx = model.nx;
    nu = model.nu;
    % x0 = [1e-1; 1e0; 2e-1; 2e0]; % pendulum
    % u = 0;

    % linear_mass_spring_
    x0 = ones(nx,1);
    u = ones(nu,1);

    %% acados sim
    sim = AcadosSim();
    sim.model.name = model_name;
    sim.model.x = model.sym_x;
    sim.model.u = model.sym_u;
    sim.model.xdot = model.sym_xdot;
    sim.solver_options.Tsim = Ts;
    sim.solver_options.integrator_type = upper(method);
    sim.solver_options.num_stages = num_stages;
    sim.solver_options.num_steps = num_steps;
    sim.solver_options.newton_iter = newton_iter;
    sim.solver_options.sens_forw = sens_forw;
    sim.solver_options.jac_reuse = jac_reuse;
    if strcmp(method, 'erk')
        sim.model.f_expl_expr = model.dyn_expr_f_expl;
    else
        sim.model.f_impl_expr = model.dyn_expr_f_impl;
    end
    sim_solver = AcadosSimSolver(sim);

    % Note: this does not work with gnsf, because it needs to be available
    % in the precomputation phase
    % 	sim_solver.set('T', Ts);

    % set initial state
    sim_solver.set('x', x0);
    sim_solver.set('u', u);

    % initialize implicit integrator
    if (strcmp(method, 'irk'))
        sim_solver.set('xdot', zeros(nx,1));
    elseif (strcmp(method, 'irk_gnsf'))
        n_out = sim_solver.sim.model.gnsf_model.dims.nout;
        sim_solver.set('phi_guess', zeros(n_out,1));
    end

    % solve
    sim_solver.solve();

    xn = sim_solver.get('xn');
    S_forw_ind = sim_solver.get('S_forw');

    %% compute forward sensitivities using finite differences
    S_forw_fd = zeros(nx, nx+nu);

    %% asymmetric finite differences
    for ii=1:nx
        dx = zeros(nx, 1);
        dx(ii) = 1.0;

        sim_solver.set('x', x0+FD_epsilon*dx);
        sim_solver.set('u', u);
        sim_solver.solve();

        xn_tmp = sim_solver.get('xn');
        S_forw_fd(:,ii) = (xn_tmp - xn) / FD_epsilon;
    end

    for ii=1:nu
        du = zeros(nu, 1);
        du(ii) = 1.0;

        sim_solver.set('x', x0);
        sim_solver.set('u', u+FD_epsilon*du);
        sim_solver.solve();

        xn_tmp = sim_solver.get('xn');
        S_forw_fd(:,nx+ii) = (xn_tmp - xn) / FD_epsilon;
    end

    %% compute & check error
    error_abs = max(max(abs(S_forw_fd - S_forw_ind)));
    disp(' ')
    disp(['integrator:  ' method]);
    disp(['error forward sensitivities (wrt finite differences):   ' num2str(error_abs)])
    disp(' ')
    if error_abs > 1e-6
        disp(['forward sensitivities error too large'])
        error(strcat('test_sens_adj FAIL: forward sensitivities error too large: \n',...
            'for integrator:\t', method));
    end
end

fprintf('\nTEST_FORWARD_SENSITIVITIES: success!\n\n');
