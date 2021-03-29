%% System Parameters
SYSTEM = 'RigidBodyMoving';
EXT_ACC = [0; 0; 0];
Q_0     = [0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0; 1];
V_0     = [1; 1; 1; 0; 1; 0; 0; 0; 0; 0; 1; 0];
MASS    = 10;
DIM     = 3;

%% Integrator
INTEGRATOR = 'GGL_VI_theta_B';
INT_PARA = [1 0.5];
DT    = 0.05;
T_0   = 0;
T_END = 10;

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