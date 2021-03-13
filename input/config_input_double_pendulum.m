%% System Parameters
SYSTEM  = 'DoublePendulum';
EXT_ACC = [0; 0; 0];
Q_0     = [1; 0; 0; 2.5; 0; 0];
V_0     = [0; 0.3; 0; 0; 0; 0];
MASS    = [1; 1.5];
DIM     = 3;

%% Integrator
INTEGRATOR = 'GGL_VI';
DT    = 0.05;
T_0   = 0;
T_END = 10;

%% Solver Method
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-09;

%% Postprocessing
shouldAnimate   = true;
plot_quantities = {'energy','energy_difference','angular_momentum','angular_momentum_difference','constraint_position','constraint_velocity'};

%% Write variables into a .mat-File
save(mfilename);