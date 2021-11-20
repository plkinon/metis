 %% METIS: A code framework for constrained mechanical systems
% author: Philipp Kinon
% date  : 02.12.2020
clearvars;
addpath(genpath(fileparts(which(mfilename))));
addpath('~/git/matlab2tikz/src');

%% METIS initialise
% Load configuration parameters from file into Metis-object
[simulation, system, integrator, solver] = Metis('config_input_heavy_top',1,1);

%% METIS solver
% Solve my System with my Solver and my Integration scheme
simulation = solver.solve(simulation, system, integrator);
%% METIS postprocessing
% Define Postprocessing from class
postprocess = Postprocess();

% Compute various postprocessing quantities
simulation = postprocess.compute(system, simulation);

% Animation of trajectory
postprocess.animation(system, simulation);

% Time-evolution of postprocessing quantites as plots
postprocess.plot(simulation);

% Export simulation results
postprocess.save(simulation);
