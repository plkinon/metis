%% System Parameters
% Name of system in /classes/System
SYSTEM = 'RedundantTwoMassSpringSystem';
% External acceleration
EXT_ACC = [0];
% Initial configuration [mass no.1; mass no.2; mass no.3, mass no.4]
Q_0 = [0; 2; 0];
% Initial velocity [mass no.1; mass no.2; mass no.3, mass no.4]
V_0 = [0.1; 0.1; 0.1];
% Masses [mass no.1; mass no.2; mass no.3, mass no.4]
MASS = [1; 1];
% Spatial dimensions
DIM = 1;

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'EML';
% Parameters of the method
INT_PARA = [NaN, NaN];
% time step size
DT = 0.1;
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
plot_quantities = {'energy', 'energy_difference','general_energy_function', 'energy_function_difference'};
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