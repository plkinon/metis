%% System Parameters
SYSTEM = 'HeavyTop';
EXT_ACC = [0; 0; 9.81];
DIM = 3;

rho = 2700;
a=1;
r=a/2;
l=3*a/4;
MASS    = rho*pi*r^2*a/3;
J1 = 3/80*MASS*(4*r^2+a^2);
J3 = 3/10*MASS*r^2;
% Initial configuration
alpha0 = pi/3;
%alpha0 = pi/2;

R0 = [1     0              0;
      0 cos(alpha0) -sin(alpha0);
      0 +sin(alpha0) cos(alpha0)];

d10 = R0*[1; 0; 0];
d20 = R0*[0; 1; 0];
d30 = R0*[0; 0; 1];
phi0 = R0*[0;0;l];
Q_0 = [phi0; d10; d20; d30];

omega_p = 10;
omega_s = MASS*EXT_ACC(3)*l/(J3*omega_p) + ((-J3+J1+MASS*l^2)/J3)*omega_p*cos(alpha0);

% Initial angular velocity
OMEGA0 = (omega_p*eye(3) + omega_s*R0)*[0;0;1];
v10 = cross(OMEGA0,d10);
v20 = cross(OMEGA0,d20);
v30 = cross(OMEGA0,d30);

V_0     = [0; 0; 0; v10; v20; v30];

clear rho a r l d J1 J3 alpha0 R0 phi0 d10 d20 d30 OMEGA0 v10 v20 v30 omega_p omega_s

%% Integrator
INTEGRATOR = 'GGL_VI_theta_B';
INT_PARA = [1 0.5];
DT    = 0.001;
T_0   = 0;
T_END = 2;

%% Solver Method
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-09;

%% Postprocessing
shouldAnimate   = true;
plot_quantities = {'energy','energy_difference','angular_momentum','angular_momentum_difference','constraint_position','constraint_velocity'};
should_export         = false;
should_export_figures = false;
export_path           = 'scratch/';

%% Write variables into a .mat-File
save(mfilename);