%% System Parameters
% Name of system in /classes/System
SYSTEM = 'Pendulum';
% External acceleration
EXT_ACC = [0; 0; -9.81];
phi0 = pi/3; l = 10;
% Initial configuration
Q_0 = [-l*sin(phi0); 0; -l*cos(phi0)];
% Initial velocity
V_0 = [0; 0; 0];
% Total mass
MASS = 1000;
% Number of spatial dimensions
DIM = 3;

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'EMS_ggl';
% Parameters of the method
INT_PARA = [NaN, NaN];
% time step size
DT = 0.01;
% starting time
T_0 = 0;
% end time
T_END = 8;

%% Solver Method
% maximum number of iterations of Newton Rhapson method
MAX_ITERATIONS = 40;
% tolerance of Newton Rhapson method
TOLERANCE = 1E-09;

clear l phi0

%% Postprocessing
% Animation of trajectory [true/false]
shouldAnimate = false;
% List of desired quantities for plotting in postprocessing
plot_quantities = {'energy', 'energy_difference', 'angular_momentum', 'angular_momentum_difference', 'constraint_position', 'constraint_velocity'};
% Export of simulation results in a .mat-file [true/false]
should_export = true;
% Export of figures in .eps- and .tikz-files
should_export_figures = true;
% Path where export-folder is created
export_path = 'scratch/';
% Matlab2Tikz (metis searches for matlab2tikz here. if not available, it
% clones the matlab2tikz repository there)
matlab2tikz_directory = '~/git/matlab2tikz';

%% Write variables into a .mat-File
% for further processing by metis
save(mfilename);