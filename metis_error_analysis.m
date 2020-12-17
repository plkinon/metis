%% METIS: A code framework for constrained mechanical systems
% author: Philipp Kinon
% date  : 02.12.2020
clearvars;
addpath(genpath(fileparts(which(mfilename))));

%% METIS initialise
% Load configuration parameters from file
CONFIG = metis_load('config_input_error_analysis');

q4 = zeros(length(CONFIG.ALL_DT),CONFIG.DIM);
error = zeros(length(CONFIG.ALL_DT)-1,1);

for i = 1:length(CONFIG.ALL_DT)
    
    CONFIG.DT = CONFIG.ALL_DT(i);

    % Initialise system, integrator and solver by class
    [system,integrator,solver] = metis_initialise(CONFIG);

    %% METIS solver
    % Solve my System with my Solver and my Integration scheme
    system = solver.solve(system,integrator);

    q4(i,:) = system.z(end,3*CONFIG.DIM+1:4*CONFIG.DIM);

end

for i = 1:length(CONFIG.ALL_DT)-1
    error(i) = abs(q4(i,:)-q4(end,:))/abs(q4(end,:));
end

figure()
loglog(CONFIG.ALL_DT(1:end-1),error)
title('EMS without GGL: relative error in q')
grid on
xlabel('h');
ylabel('relative error');
