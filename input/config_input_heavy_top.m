%% System Parameters
SYSTEM = 'HeavyTop';
EXT_ACC = [0; 0; 9.81];
DIM = 3;

rho = 2700;
a=0.1;
r=a/2;
l=3*a/4;
d=l;
MASS    = rho*pi*r^2*a/3;
J2 = 3/80*MASS*(4*r^2+a^2) + MASS*d^2;
J3 = 3/10*MASS*r^2;
% Initial configuration
alpha0 = pi/3;

R0 = [1 0 0;
      0 cos(alpha0) -sin(alpha0);
      0 sin(alpha0) cos(alpha0)];

d10 = R0*[1; 0; 0];
d20 = R0*[0; 1; 0];
d30 = R0*[0; 0; 1];
Q_0     = [0; 0; 0; d10; d20; d30];

phi0 = 10;
psi0 = MASS*EXT_ACC(3)*l/(J3*phi0) + (J2-J3)/J3*phi0*cos(alpha0);
% Initial angular velocity
OMEGA0 = [0; phi0*sin(alpha0); psi0+phi0*cos(alpha0)];
v10 = cross(OMEGA0,d10);
v20 = cross(OMEGA0,d20);
v30 = cross(OMEGA0,d30);
V_0     = [0; 0; 0; v10; v20; v30];

clear rho a r l d J2 J3 alpha0 R0 d10 d20 d30 phi0 psi0 OMEGA0 v10 v20 v30
%% Integrator
INTEGRATOR = 'GGL_VI_theta_B';
INT_PARA = [1 0.5];
DT    = 0.001;
T_0   = 0;
T_END = 1;

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