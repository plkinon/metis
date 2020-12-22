%% Problem Parameters
SYSTEM = 'Pendulum';
EXT_ACC = [0; 0; 9.81];
Q_0     = [1; 0; 0];
V_0     = [0; 0.2; 0.1];
MASS    = 1;
DIM     = 3;
nDOF    = 3;

%% Integrator
INTEGRATOR = 'Ggl_variational';
DT    = 0.02;
T_0   = 0;
T_END = 5;

%% Solver Method
SOLVER         = 'Newton';
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-09;

%% Postprocessing
shouldAnimate   = true;
plot_quantities = {'energy','energy_difference','angular_momentum','angular_momentum_difference','constraint_position','constraint_velocity'};

%% Write variables into a .mat-File
save(mfilename);