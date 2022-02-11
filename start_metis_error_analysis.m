%% metis_error_analysis.m - starts error analysis with metis
%
% Metis is an object-oriented MATLAB code package for the simulation of
% constrained mechanical systems under the usage of numerical
% time-integration methods and Newton-Rhapson method.
%
% Usage:
%       metis_error_analysis()
%
% Other .m-files required: error-analysis-input-file in /input
% .mat-files required: none
%
% Author : Philipp Kinon
% Email  : philipp.kinon@kit.edu
% Date   : 02.12.2021

%% ----------------------------BEGIN CODE ---------------------------------

%% METIS initialise
% Clear present variables
clearvars;
% Add all subdirectories and matlab2tikz to the current path
addpath(genpath(fileparts(which(mfilename))));

% Metis creates a dummy simulation object, the system and solver from input-file
[dummy_simulation, system, ~, solver] = Metis('input/published/paper_kinon_betsch_schneider_2022/error_analysis_heavy_top', 1, 1);
% Check how many different timestepsizes and integrators are analyzed
n_DT = numel(dummy_simulation.ALL_DT);
n_INT = numel(dummy_simulation.ALL_INTEGRATOR);

% Define Postprocessing from class
postprocess = Postprocess();
dummy_simulation = postprocess.compute(system, dummy_simulation);

% Analysis quantity allocation
analyzed_quantity = zeros(system.DIM, n_DT, n_INT);

%% Loop over all desired timestepsizes and integration schemes
for i = 1:n_DT
    for j = 1:n_INT

        % Metis creates objects for current timestepsize and integrator
        [current_simulation, ~, current_integrator, ~] = Metis('input/published/paper_kinon_betsch_schneider_2022/error_analysis_heavy_top', i, j);

        %% METIS solver
        % Solve system with solver and current integrator
        current_simulation = solver.solve(current_simulation, system, current_integrator);

        % Compute various postprocessing quantities
        current_simulation = postprocess.compute(system, current_simulation);

        % Compute position error and velocity constraint violation
        analyzed_quantity(:, i, j) = system.hconvergence_set(current_simulation);

    end
end

%% Error computation
% Compute relative error w.r.t. to smalles h (last entry in ALL_DT)
reference_solution = system.hconvergence_reference(dummy_simulation, analyzed_quantity);
error              = postprocess.calculate_errors(analyzed_quantity,reference_solution,n_DT,n_INT);

%% Plot and export
postprocess.convergence_plot(dummy_simulation.ALL_DT, error, n_INT, dummy_simulation.ALL_INTEGRATOR);
addpath([current_simulation.matlab2tikz_directory,'/src']);
matlab2tikz('height', '\figH', 'width', '\figW', 'filename', 'scratch/hconvergence.tikz', 'showInfo', false, 'floatformat', '%.5g');

% -------------------------- END OF CODE ----------------------------------