%% METIS: A code framework for constrained mechanical systems
% author: Philipp Kinon
% date  : 02.05.2021
clearvars;
addpath(genpath(fileparts(which(mfilename))));
addpath('~/git/matlab2tikz/src');

%% METIS initialise
% Load configuration parameters from file into Metis-object
[dummy_simulation,system,~,solver] = Metis('config_input_error_analysis_heavy_top',1,1);
n_DT = numel(dummy_simulation.ALL_DT);
n_INT = numel(dummy_simulation.ALL_INTEGRATOR);

% Define Postprocessing from class
postprocess = Postprocess();
dummy_simulation = postprocess.compute(system,dummy_simulation);

% Analysis output quantity allocation
analyzed_quantity  = zeros(system.DIM,n_DT,n_INT);

%% Loop over all desired timestepsizes and integration schemes
for i = 1:n_DT
    for j = 1:n_INT
    
    % Load configuration parameters for current timestepsize
    [current_simulation,~,current_integrator,~] = Metis('config_input_error_analysis_heavy_top',i,j);
     
    %% METIS solver
    % Solve my System with my Solver and my Integration scheme
    current_simulation = solver.solve(current_simulation,system,current_integrator);

    % Compute various postprocessing quantities
    current_simulation = postprocess.compute(system,current_simulation);
    
    % Compute position error and velocity constraint violation
    analyzed_quantity(:,i,j) = system.hconvergence_set(current_simulation);
        
    end
end

%% Compute relative error w.r.t. to smalles h (last entry in ALL_DT)
reference_solution = system.hconvergence_reference(dummy_simulation,analyzed_quantity);
error              = postprocess.calculate_errors(analyzed_quantity,reference_solution,n_DT,n_INT);

%% Figure: positional error
postprocess.convergence_plot(dummy_simulation.ALL_DT,error,n_INT,dummy_simulation.ALL_INTEGRATOR);
matlab2tikz('height','\figH','width','\figW','filename','scratch/hconvergence.tikz','showInfo', false,'floatformat','%.5g');
