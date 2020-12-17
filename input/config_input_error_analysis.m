%% Problem Parameters
SYSTEM  = 'FourParticleSystem';
EXT_ACC = [0; 0; 0];
Q_0     = [0; 0; 0; 1; 0; 0; 0; 1; 0; 1; 1; 0];
V_0     = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2/1.7];
MASS    = [1; 3; 2.3; 1.7];
DIM     = 3;

%% Integrator
INTEGRATOR = 'EMS_std';
INTEGRATION_VARIABLES = 0.5;
ALL_DT    = [1 0.1 0.01 0.001];
T_0   = 0;
T_END = 10;

%% Solver Method
SOLVER         = 'Newton';
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-08;
NUM_TANGENT    = true;

%% Postprocessing
% no postprocessing

%% Write variables into a .mat-File
save(mfilename);