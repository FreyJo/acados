%% Simulink example
%

%% Run minimal example
%
minimal_example_ocp;


%% get available simulink_opts with default options
simulink_opts = get_acados_simulink_opts;

% manipulate simulink_opts

% inputs
simulink_opts.inputs.cost_W_0 = 1;
simulink_opts.inputs.cost_W = 1;
simulink_opts.inputs.cost_W_e = 1;

% outputs
% simulink_opts.outputs.utraj = 1;
% simulink_opts.outputs.xtraj = 1;

simulink_opts.samplingtime = '-1';
    % 't0' (default) - use time step between shooting node 0 and 1
    % '-1' - inherit sampling time from other parts of simulink model


%% Render templated Code for the model contained in ocp object
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

%%
disp('Press play in Simulink!');
