function [this_system, this_integrator, this_solver] = metis_initialise(CONFIG)
%% Function: Initialises problem, integrator and solver objects
% Author: Philipp Kinon
% date: 03.12.2020

% Clear workspace, close all windows and clear command window 
close all; clc;

fprintf('************************************************** \n');
fprintf(' METIS - Computing constrained mechanical systems \n');
fprintf('************************************************** \n');

%% Check if user input is valid
check_user_input(CONFIG);

%% System
% Takes user-defined string and evaluates it as the constructor of a class
% system
this_system = feval(CONFIG.SYSTEM,CONFIG);

%% Integrator
% Define Integrator from Class (same procedure as for system)
this_integrator = feval(CONFIG.INTEGRATOR,CONFIG);

%% Apply integrator to system and vice versa
% Initialise integrator for current system
this_integrator = this_integrator.initialise(CONFIG,this_system);

% Intialise System for current integrator
this_system = this_system.initialise(CONFIG,this_integrator);

%% Solver
% Define Solver from class
this_solver = feval(CONFIG.SOLVER,CONFIG);

end