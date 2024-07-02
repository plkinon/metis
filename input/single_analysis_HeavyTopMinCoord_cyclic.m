%% System Parameters
% Name of system in /classes/System
SYSTEM = 'HeavyTopMinCoordCyclic';
% External acceleration
EXT_ACC = 9.81;
g= EXT_ACC;
% Initial configuration
Q_0 = [pi/3; 0;0];

rho = 2700; % mass density
a = 0.1; % length of the gyro top
r = a / 2; % radius of the gyro top
l = 3 * a / 4; % location of center of mass along sym. axis
MASS = rho * pi * r^2 * a / 3; % total mass of the gyro top
J1 = 3 / 80 * MASS * (4 * r^2 + a^2); % inertia moment w.r.t. d1-axis (J1 = J2)
J3 = 3 / 10 * MASS * r^2; % inertia moment w.r.t. d3-axis (sym. axis)

%% Initial angular velocity
% precession rate (angular velocity about e3-axis)
omega_precession = 10;
alpha0 = Q_0(1);
% spin rate (angular velocity about d3-axis) for steady precession
omega_spin = MASS * g * l / (J3 * omega_precession) + ((J1 + MASS * l^2 - J3) / J3) * omega_precession * cos(alpha0);
% Initial velocity
V_0 = [0; omega_spin; omega_precession];

clear J1 J3 a r l omega_precession omega_spin alpha0 g rho
% Spatial dimensions
DIM = 3;

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'MP_noCons'; 
%INTEGRATOR = 'EMG_noCons'; % DG due to Gonzalez without accounting for cyclic coordinate
%INTEGRATOR = 'EMS_std_cyclic'; %method by DS, PK

% Parameters of the method
INT_PARA = [NaN, NaN];
% time step size
DT = 0.002;
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
shouldAnimate = true;
% List of desired quantities for plotting in postprocessing
plot_quantities = {'energy', 'energy_difference','general_energy_function', 'energy_function_difference', 'angular_momentum', 'angular_momentum_difference'};
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