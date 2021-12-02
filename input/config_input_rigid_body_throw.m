%% System Parameters
% Name of system in /classes/System
SYSTEM = 'RigidBodyMoving';
% External acceleration
EXT_ACC = [0; 0; -10];
% Initial configuration [center of mass; director no.1; no.2; no.3]
Q_0     = [0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0; 1];
% Initial velocity [center of mass; director no.1; no.2; no.3]
V_0     = [10; 10; 10; 0; 1; -1; -1; 0; 0; 1; 0; 0];
% Total mass
MASS    = 1;
% Number of spatial dimensions
DIM     = 3;

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'EMS_ggl';
% Parameters of the method
INT_PARA = [NaN NaN];
% time step size
DT    = 0.02;
% starting time
T_0   = 0;
% end time
T_END = 2;

%% Solver Method
% maximum number of iterations of Newton Rhapson method
MAX_ITERATIONS = 40;
% tolerance of Newton Rhapson method
TOLERANCE      = 1E-09;

%% Postprocessing
% Animation of trajectory [true/false]
shouldAnimate = false;
% List of desired quantities for plotting in postprocessing
plot_quantities = {'energy','energy_difference','angular_momentum','angular_momentum_difference','constraint_position','constraint_velocity'};
% Export of simulation results in a .mat-file [true/false]
should_export = false;
% Export of figures in .eps- and .tikz-files
should_export_figures = false;
% Path where export-folder is created
export_path = 'scratch/';

%% Write variables into a .mat-File
% for further processing by metis
save(mfilename);