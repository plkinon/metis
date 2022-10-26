%% System Parameters
% Name of system in /classes/System
SYSTEM = 'FourParticleSystem';
% External acceleration
EXT_ACC = [0; 0; 0];
% Initial configuration [mass no.1; mass no.2; mass no.3, mass no.4]
Q_0 = [0; 0; 0; 1; 0; 0; 0; 1; 0; 1; 1; 0];
% Initial velocity [mass no.1; mass no.2; mass no.3, mass no.4]
V_0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 3 / 1.7; 2 / 1.7];
% Masses [mass no.1; mass no.2; mass no.3, mass no.4]
MASS = [1; 3; 2.3; 1.7];
% Spatial dimensions
DIM = 3;

%% Integrator
% Name of routines in /classes/Integrator to be analyzed
INTEGRATOR = {'EMS_ggl', 'EMS_std'};
% corresponding parameters
INT_PARA = [NaN, NaN; NaN, NaN];
% time step sizes to be analyzed
DT = [0.1, 0.01, 0.001, 0.0001];
% quantity which is to be analyzed
CONV_QUANTITY = 'p';
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