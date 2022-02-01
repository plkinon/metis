%% System Parameters
% Name of system in /classes/System
SYSTEM = 'Pendulum';
% External acceleration
EXT_ACC = [0; 0; 9.81];
% Initial configuration
Q_0 = [1; 0; 0];
% Initial velocity
V_0 = [0; 1; 0];
% Total mass
MASS = 1;
% Number of spatial dimensions
DIM = 3;

%% Integrator
% Name of routines in /classes/Integrator to be analyzed
INTEGRATOR = {'CSE_B', 'GGL_std', 'GGL_VI'};
% corresponding parameters
INT_PARA = [NaN, NaN; NaN, NaN; NaN, NaN];
% time step sizes to be analyzed
DT = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001];
% quantity which is to be analyzed
CONV_QUANTITY = 'q';
% starting time
T_0 = 0;
% end time
T_END = 1;

%% Solver Method
% maximum number of iterations of Newton Rhapson method
MAX_ITERATIONS = 40;
% tolerance of Newton Rhapson method
TOLERANCE = 1E-09;

%% Postprocessing
% Animation of trajectory [true/false]
shouldAnimate = false;
% List of desired quantities for plotting in postprocessing
plot_quantities = {};
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