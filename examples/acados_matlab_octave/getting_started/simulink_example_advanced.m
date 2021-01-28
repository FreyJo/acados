%% Simulink example
%

%% Run minimal example
%
minimal_example_ocp;


%% Render templated Code for the model contained in ocp object
%
simulink_opts = get_acados_simulink_opts;

% manipulate simulink_opts to [de]activate in- & outputs (preliminary)
simulink_opts.inputs.cost_W_0 = 1;
simulink_opts.inputs.cost_W = 1;
simulink_opts.inputs.cost_W_e = 1;

% simulink_opts.outputs.sqp_iter = 0;
% simulink_opts.outputs.utraj = 1;
% simulink_opts.outputs.xtraj = 1;



ocp.generate_c_code(simulink_opts);

%% Compile Sfunctions
cd c_generated_code

make_sfun_sim; % integrator
make_sfun; % ocp solver


%% Copy Simulink example blocks into c_generated_code
source_folder = fullfile(pwd, '..');
target_folder = pwd;
copyfile( fullfile(source_folder, 'simulink_model_advanced_closed_loop.slx'), target_folder );

%% Open Simulink example blocks
open_system(fullfile(target_folder, 'simulink_model_advanced_closed_loop'))

% copyfile( fullfile(source_folder, 'simulink_model_integrator.slx'), target_folder );
% open_system(fullfile(target_folder, 'simulink_model_integrator'))


%%
disp('Press play in Simulink!');
