%% Problem Parameters
SYSTEM = 'Pendulum';
EXT_ACC = [0; 0; 9.81];
Q_0     = [1; 0; 0];
V_0     = [0; 0.2; 0];
DIM     = 3;
nDOF    = 3;

%% Integrator
INTEGRATOR = 'Ggl_std';
INTEGRATION_VARIABLES = 0.5;
DT    = 0.02;
T_0   = 0;
T_END = 5;

%% Solver Method
SOLVER         = 'Newton';
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-08;
NUM_TANGENT    = false;

%% Postprocessing
shouldAnimate   = true;
plot_quantities =  {'energy','constraint_position'};

%% Write variables into a .mat-File
save(mfilename);