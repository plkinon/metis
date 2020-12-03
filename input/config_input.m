%% Time-stepping
DT    = 0.02;
T_0   = 0;
T_END = 5;

%% Integrator
INTEGRATION_VARIABLES = 0.5;

%% Newton-Rhapson Method
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-08;
NUM_TANGENT    = false;

%% Physical Parameters
EXT_ACC = [0; 0; 9.81];
Q_0     = [1; 0; 0];
V_0     = [0; 0.2; 0];
DIM     = 3;
nDOF    = 3;

%% Postprocessing
shouldAnimate   = true;
plot_quantities =  {'energy','constraint_position'};

%% Write variables into a .mat-File
save(mfilename);