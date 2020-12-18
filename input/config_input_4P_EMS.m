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
DT    = 0.02;
T_0   = 0;
T_END = 10;

%% Solver Method
SOLVER         = 'Newton';
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-08;
NUM_TANGENT    = true;

%% Postprocessing
shouldAnimate   = true;
plot_quantities = {'energy','energy_difference','angular_momentum','angular_momentum_difference','constraint_position','constraint_velocity'};

%% Write variables into a .mat-File
save(mfilename);