%% METIS: A code framework for constrained mechanical systems
% author: Philipp Kinon
% date  : 02.05.2021
clearvars;
addpath(genpath(fileparts(which(mfilename))));
addpath('~/git/matlab2tikz/src');

%% METIS initialise
% Load configuration parameters from file into Metis-object
[dummy_simulation,system,~,solver] = Metis('config_input_error_analysis',1,1);
n_DT = numel(dummy_simulation.ALL_DT);
n_INT = numel(dummy_simulation.ALL_INTEGRATOR);

% Define Postprocessing from class
postprocess = Postprocess();
dummy_simulation = postprocess.compute(system,dummy_simulation);

% Analysis output quantities
qEnd  = zeros(size(dummy_simulation.Q_0,1),numel(dummy_simulation.ALL_DT),numel(dummy_simulation.ALL_INTEGRATOR));   %position 

%% Loop over all desired timestepsizes and integration schemes
for i = 1:length(dummy_simulation.ALL_DT)
    for j = 1:length(dummy_simulation.ALL_INTEGRATOR)
    
    % Load configuration parameters for current timestepsize
    [current_simulation,~,current_integrator,~] = Metis('config_input_error_analysis',i,j);
     
    %% METIS solver
    % Solve my System with my Solver and my Integration scheme
    current_simulation = solver.solve(current_simulation,system,current_integrator);

    % Compute various postprocessing quantities
    current_simulation = postprocess.compute(system,current_simulation);
    
    % Compute position error and velocity constraint violation
    qEnd(:,i,j) = current_simulation.z(end,1:system.nBODIES*system.DIM);
        
    end
end

%% Compute relative error w.r.t. to smalles h (last entry in ALL_DT)
error = postprocess.calculate_errors(qEnd,n_DT,n_INT,n_DT);

%% Figure: positional error
postprocess.convergence_plot(dummy_simulation.ALL_DT,error,n_INT,dummy_simulation.ALL_INTEGRATOR);