%% System Parameters
SYSTEM = 'HeavyTop';
g = 9.81;
EXT_ACC = [0; 0; -g];
DIM = 3;

%% Geometric parameters
rho  = 2700;                    % mass density
a    = 0.1;                       % length of the gyro top
r    = a/2;                     % radius of the gyro top
l    = 3*a/4;                   % location of center of mass along sym. axis
MASS = rho*pi*r^2*a/3;          % total mass of the gyro top
J1   = 3/80*MASS*(4*r^2+a^2);   % inertia moment w.r.t. d1-axis (J1 = J2)
J3   = 3/10*MASS*r^2;           % inertia moment w.r.t. d3-axis (sym. axis)

%% Initial configuration
% Initial inclination within e2-e3-plane
alpha0 = pi/3;

% Rotation matrix for rotation about e1-axis
R0 = [1     0              0;
      0 cos(alpha0) -sin(alpha0);
      0 +sin(alpha0) cos(alpha0)];

% Transformation of directors and center of mass
d10 = R0*[1; 0; 0];
d20 = R0*[0; 1; 0];
d30 = R0*[0; 0; 1];
phi0 = R0*[0;0;l];
Q_0 = [phi0; d10; d20; d30];


%% Initial angular velocity
% precession rate (angular velocity about e3-axis)
omega_p = 10;
% spin rate (angular velocity about d3-axis) for steady precession
omega_s = MASS*g*l/(J3*omega_p) + ((J1+MASS*l^2-J3)/J3)*omega_p*cos(alpha0);
% inital angular velocity vector w.r.t. e_i-coordinate system
OMEGA0 = (omega_p*eye(3) + omega_s*R0)*[0;0;1];

%% Initial velocity vector
% initial velocity of center of mass
v00 = cross(OMEGA0,phi0);
% initial velocities of directors
v10 = cross(OMEGA0,d10);
v20 = cross(OMEGA0,d20);
v30 = cross(OMEGA0,d30);
V_0     = [v00; v10; v20; v30];

% clear unnecessary variables
clear rho a r l d g J1 J3 alpha0 R0 phi0 d10 d20 d30 OMEGA0 v00 v10 v20 v30 omega_p omega_s

%% Integrator
INTEGRATOR = {'GGL_VI_mod'};
INT_PARA = [NaN NaN];
DT = [0.0001 0.00005 0.00001 0.000005];
T_0   = 0;
T_END = 0.001;

%% Solver Method
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-09;

%% Postprocessing
shouldAnimate   = false;
plot_quantities = {'energy','energy_difference','angular_momentum','angular_momentum_difference','constraint_position','constraint_velocity'};
should_export         = false;
should_export_figures = false;
export_path           = 'scratch/';

%% Write variables into a .mat-File
save(mfilename);