%% METIS: A code framework for constrained mechanical systems
% author: Philipp Kinon
% date  : 02.12.2020
clearvars;
addpath(genpath(fileparts(which(mfilename))));

%% METIS initialise
% Load configuration parameters from file into Metis-object
[simulation,system,integrator,solver] = Metis('config_input_error_analysis');

% Define Postprocessing from class
postprocess = Postprocess();

% Analysis output quantities
qEnd  = zeros(length(simulation.ALL_DT),size(simulation.Q_0,1));   %position 
error = zeros(length(simulation.ALL_DT),4);                        %position error w.r.t. to smallest h
consV = zeros(length(simulation.ALL_DT),1);                        %velocity constraint violation

%% Loop over all desired timestepsizes
for i = 1:length(simulation.ALL_DT)
    
    % Apply current timestepsize
    simulation.DT = simulation.ALL_DT(i);
    integrator.DT = simulation.ALL_DT(i);
    integrator.t  = simulation.T_0:simulation.DT:simulation.T_END;
    simulation    = simulation.set_solution_matrix(integrator,system);
    
    %% METIS solver
    % Solve my System with my Solver and my Integration scheme
    simulation = solver.solve(simulation,system,integrator);

    % Compute various postprocessing quantities
    simulation = postprocess.compute(system,simulation);
    
    % Compute position error and velocity constraint violation
    qEnd(i,:) = simulation.z(end,1:system.nBODIES*system.DIM);
    consV(i)  = norm(simulation.constraint_velocity(end));
    
end

%% Compute relative error w.r.t. to smalles h
for i = 1:length(simulation.ALL_DT)
    for j = 1:system.nBODIES
        DIM = system.DIM;
        error(i,j) = abs(qEnd(i,(j-1)*DIM+1:j*DIM)-qEnd(end,(j-1)*DIM+1:j*DIM))/abs(qEnd(end,(j-1)*DIM+1:j*DIM));
    end
end

%% Figure 1: position error
figure()
loglog(simulation.ALL_DT(1:end-1),error(1:end-1,:),'Linewidth',2)
title([integrator.NAME,': relative error in q'])
grid on
xlabel('h');
ylabel('relative error');

%% Figure 2: velocity constraint violation (norm)
figure()
loglog(simulation.ALL_DT,consV,'Linewidth',2)
title([integrator.NAME,' : velocity constraint norm'])
grid on
xlabel('h');
ylabel('norm(G*v)');