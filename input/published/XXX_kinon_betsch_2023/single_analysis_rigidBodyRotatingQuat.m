%% System Parameters
% Name of system in /classes/System
SYSTEM = 'RigidBodyRotatingQuaternions';
%SYSTEM = 'RigidBodyRotatingQuaternionsRegularMassMatrix';

% External acceleration
EXT_ACC = 0;
% Initial configuration
Q_0 = [1; 0; 0; 0];
% Initial velocity 
Omega_0 = [10;20;20];
%extract vector and scalar part form quaternion
Q0_vec = Q_0(2:4);
Q0_scalar = Q_0(1);

%skew-sym matrix corresponding to vector part
Q0_hat = [0, -Q0_vec(3), Q0_vec(2);
        Q0_vec(3), 0, -Q0_vec(1);
        -Q0_vec(2), Q0_vec(1), 0];

% transformation matrix
G_Q0 = [-Q0_vec, Q0_scalar*eye(3) - Q0_hat];
V_0 = 1/2*G_Q0'*Omega_0;
% Mass
MASS = 1;
% Spatial dimensions
DIM = 3;

% clear unnecessary variables (crucial for further processing!)
clear Omega_0 Q0_hat Q0_vec Q0_scalar G_Q0

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'EML';
%INTEGRATOR = 'EMS_std';

% Parameters of the method
INT_PARA = [NaN, NaN];
% time step size
DT = 0.05;
% starting time
T_0 = 0;
% end time
T_END = 2;

%% Solver Method
% maximum number of iterations of Newton Rhapson method
MAX_ITERATIONS = 40;
% tolerance of Newton Rhapson method
TOLERANCE = 1E-09;

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