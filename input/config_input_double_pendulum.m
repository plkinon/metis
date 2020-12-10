%% Problem Parameters
SYSTEM  = 'DoublePendulum';
EXT_ACC = [0; 0; 9.81];
Q_0     = [1; 0; 0; 2.5; 0; 0];
V_0     = [0; 0.3; 0; 0; 0; 0];
MASS    = [1; 1.5];
DIM     = 3;

%% Integrator
INTEGRATOR = 'Ggl_gen_alpha';
INTEGRATION_VARIABLES = 0.5;
DT    = 0.02;
T_0   = 0;
T_END = 10;

%% Solver Method
SOLVER         = 'Newton';
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-08;
NUM_TANGENT    = false;

%% Postprocessing
shouldAnimate   = false;
plot_quantities =  {'energy'};

%% Write variables into a .mat-File
save(mfilename);