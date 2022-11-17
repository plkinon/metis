%% System Parameters
% Name of system in /classes/System
SYSTEM = 'DoubleFourBarLinkage';
% External acceleration
EXT_ACC = [0; -9.81];
% Initial configuration [bar no.1; bar no.2; bar no.3; bar no.4; bar no.5]
l = 1;
q10 = [0;    l/2; 0;  1;  1;  0];
q20 = [l/2;    l; 1;  0;  0; -1];
q30 = [l;    l/2; 0; -1; -1;  0];
q40 = [3*l/2;  l; 1;  0;  0; -1];
q50 = [2*l;  l/2; 0; -1; -1;  0];
Q_0 = [q10; q20; q30; q40; q50]; 
% Initial velocity [bar no.1; bar no.2; bar no.3; bar no.4; bar no.5]
v0 = 1;
v10 = [v0/2; 0; 0; 0; 0; 0];
v20 = [v0;   0; 0; 0; 0; 0];
v30 = [v0/2; 0; 0; 0; 0; 0];
v40 = [v0;   0; 0; 0; 0; 0];
v50 = [v0/2; 0; 0; 0; 0; 0];
V_0 = [v10; v20; v30; v40; v50]; 
% Masses [mass no.1; mass no.2; mass no.3, mass no.4]
MASS = [1; 1; 1; 1; 1];
% Spatial dimensions
DIM = 2;

% clear unnecessary variables (crucial for further processing!)
clear l q10 q20 q30 q40 q50 v0 v10 v20 v30 v40 v50

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'GGL_VI_mod';
% Parameters of the method
INT_PARA = [NaN, NaN];
% time step size
DT = 0.001;
% starting time
T_0 = 0;
% end time
T_END = 10;

%% Solver Method
% maximum number of iterations of Newton Rhapson method
MAX_ITERATIONS = 40;
% tolerance of Newton Rhapson method
TOLERANCE = 1E-09;

%% Postprocessing
% Animation of trajectory [true/false]
shouldAnimate = false;
% List of desired quantities for plotting in postprocessing
plot_quantities = {'energy', 'energy_difference', 'general_energy_function', 'energy_function_difference', 'angular_momentum', 'angular_momentum_difference', 'constraint_position', 'constraint_velocity'};
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