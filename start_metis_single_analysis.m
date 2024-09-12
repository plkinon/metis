%% metis_start.m - starts single simulations with metis
%
% Metis is an object-oriented MATLAB code package for the simulation of
% constrained mechanical systems under the usage of numerical
% time-integration methods and Newton-Rhapson method.
%
% This startscript provides single simulations. For error analyses see:
% metis_error_analysis.m
%
% Usage :
%       metis_start()
%
% Other .m-files required: input-file in /input
% .mat-files required: none
%
% Author : Philipp Kinon
% Email  : philipp.kinon@kit.edu
% Date   : 04.11.2022

%% ----------------------------BEGIN CODE ---------------------------------

%% METIS initialise
% Clear present variables
clearvars;
% Add all subdirectories and to the current path
addpath(genpath(fileparts(which(mfilename))));
% Metis creates objects from input-file
[simulation, system, integrator, solver] = Metis('input/single_analysis_4P', 1, 1);

%% METIS solver
% Solve system with chosen solver and integration scheme
simulation = solver.solve(simulation, system, integrator);

%% METIS postprocessing
% Define postprocessing from class
postprocess = Postprocess();

% Compute various postprocessing quantities
simulation = postprocess.compute(system, simulation, integrator);

% Animation of trajectory if activated in input-file
postprocess.animation(system, simulation);

% Plot time-evolution of postprocessing quantites
postprocess.plot(simulation);

% Export simulation results if activated in input-file
postprocess.save(simulation);

% -------------------------- END OF CODE ----------------------------------