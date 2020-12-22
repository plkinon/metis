%% METIS: A code framework for constrained mechanical systems
% author: Philipp Kinon
% date  : 02.12.2020
clearvars;
addpath(genpath(fileparts(which(mfilename))));

%% METIS initialise
% Load configuration parameters from file into Metis-object
simulation = Metis('config_input_4P_EMS');

% Initialise system, integrator and solver by class
[system,integrator,solver] = simulation.initialise();

%% METIS solver
% Solve my System with my Solver and my Integration scheme
system = solver.solve(system,integrator);

%% METIS postprocessing
% Define Postprocessing from class
postprocess = Postprocess();

% Compute various postprocessing quantities
system = postprocess.compute(system);

% Animation of trajectory
postprocess.animation(system,integrator,simulation);

% Time-evolution of postprocessing quantites as plots
postprocess.plot(simulation,system,integrator);