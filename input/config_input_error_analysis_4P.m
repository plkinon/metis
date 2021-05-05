%% System Parameters
SYSTEM  = 'FourParticleSystem';
EXT_ACC = [0; 0; 0];
Q_0     = [0; 0; 0; 1; 0; 0; 0; 1; 0; 1; 1; 0];
V_0     = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 3/1.7; 2/1.7];
MASS    = [1; 3; 2.3; 1.7];
DIM     = 3;

%% Integrator
INTEGRATOR = {'EMS_std','EMS_ggl'};
DT    = [0.1 0.01 0.001 0.0001];
INT_PARA = [NaN NaN; NaN NaN];
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
export_path           = 'scratch/EMS';

%% Write variables into a .mat-File
save(mfilename);