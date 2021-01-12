%% System Parameters
SYSTEM  = 'FourParticleSystem';
EXT_ACC = [0; 0; 0];
Q_0     = [0; 0; 0; 1; 0; 0; 0; 1; 0; 1; 1; 0];
V_0     = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 20/1.7];
MASS    = [1; 3; 2.3; 1.7];
DIM     = 3;

%% Integrator
INTEGRATOR = 'EMS_ggl';
ALL_DT = [0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002 0.0001];
T_0   = 0;
T_END = 1;

%% Solver Method
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-09;

%% Postprocessing
% no postprocessing

%% Write variables into a .mat-File
save(mfilename);