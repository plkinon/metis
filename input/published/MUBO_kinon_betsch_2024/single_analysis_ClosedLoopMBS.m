%% System Parameters
% Name of system in /classes/System
SYSTEM = 'ClosedLoopMBSQuaternions';

% External acceleration
EXT_ACC = 0;
% Initial configuration
length_bar = 10;
phi_1 = [length_bar/2; 0; 0];
phi_2 = [0; length_bar/2; 0];
phi_3 = [-length_bar/2; 0; 0];
phi_4 = [0; -length_bar/2; 0];
q_1 = [1; 0; 0; 0];
q_2 = [1; 0; 0; 0];
q_3 = [1; 0; 0; 0];
q_4 = [1; 0; 0; 0];

Q_0 = [phi_1; q_1; phi_2; q_2; phi_3; q_3; phi_4; q_4];
% Initial velocity 
V_0 = zeros(size(Q_0));
% Mass
rho = 1;
cross_section = 1;
MASS = rho * length_bar * cross_section; 
% Spatial dimensions
DIM = 3;

% clear unnecessary variables (crucial for further processing!)
clear length_bar phi_1 phi_2 phi_3 phi_4 q_1 q_2 q_3 q_4 rho cross_section

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'EML';
%INTEGRATOR = 'MP_Livens';

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
TOLERANCE = 1E-9;

%% Postprocessing
% Animation of trajectory [true/false]
shouldAnimate = false;
% List of desired quantities for plotting in postprocessing
plot_quantities = {'energy', 'energy_difference','general_energy_function', 'energy_function_difference', 'angular_momentum', 'angular_momentum_difference', 'constraint_velocity', 'constraint_position'};
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