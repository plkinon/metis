%% System Parameters
% Name of system in /classes/System
SYSTEM = 'HeavyTopQuaternions';
%SYSTEM = 'HeavyTopQuaternionsRegularMassMatrix';

% External acceleration
g = 9.81;
EXT_ACC = [0; 0; -g];
%% Geometric parameters
rho = 2700; % mass density
H = 0.1; % length of the gyro top
R = 0.05; % radius of the gyro top
L = 3 * H / 4; % location of center of mass along sym. axis
MASS = rho * pi * R^2 * H / 3; % total mass of the gyro top
J1 = 3 / 80 * MASS * (4 * R^2 + H^2); % inertia moment w.r.t. d1-axis (J1 = J2)
J3 = 3 / 10 * MASS * R^2; % inertia moment w.r.t. d3-axis (sym. axis)

% Initial configuration
theta_0 = pi/3;
Q_0 = [cos(theta_0/2); sin(theta_0/2); 0; 0];
% Rotation matrix for rotation about e1-axis
R0 = [1, 0, 0; 0, cos(theta_0), -sin(theta_0); 0, +sin(theta_0), cos(theta_0)];
% Initial velocity 
% precession rate (angular velocity about e3-axis)
omega_p = 10;
% spin rate (angular velocity about d3-axis) for steady precession
omega_s = MASS * g * L / (J3 * omega_p) + ((J1 + MASS * L^2 - J3) / J3) * omega_p * cos(theta_0);
% inital angular velocity vector w.r.t. e_i-coordinate system
omega_0 = (omega_p * eye(3) + omega_s * R0) * [0; 0; 1];
%extract vector and scalar part form quaternion
Q0_vec = Q_0(2:4);
Q0_scalar = Q_0(1);

%skew-sym matrix corresponding to vector part
Q0_hat = [0, -Q0_vec(3), Q0_vec(2);
        Q0_vec(3), 0, -Q0_vec(1);
        -Q0_vec(2), Q0_vec(1), 0];

% transformation matrix
E_Q0 = [-Q0_vec, Q0_scalar*eye(3) + Q0_hat];
V_0 = 1/2*E_Q0'*omega_0;
% Spatial dimensions
DIM = 3;

% clear unnecessary variables (crucial for further processing!)
clear omega_0 Q0_hat Q0_vec Q0_scalar E_Q0 omega_s omega_p theta_0 J1 J3 L R H rho g R0

%% Integrator
% Name of routine in /classes/Integrator
INTEGRATOR = 'EML';
%INTEGRATOR = 'MP_Livens';
%INTEGRATOR = 'EMS_std';
% Parameters of the method
INT_PARA = [NaN, NaN];
% time step size
DT = 0.02;
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
plot_quantities = {'energy', 'energy_difference','general_energy_function', 'energy_function_difference', 'angular_momentum', 'angular_momentum_difference', 'constraint_velocity', 'constraint_position','cartesian_coordinates_center_of_mass'};
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