%% System Parameters
SYSTEM  = 'Pendulum';
EXT_ACC = [0; 0; 9.81];
Q_0     = [1; 0; 0];
V_0     = [0; 1; 0];
MASS    = 1;
DIM     = 3;

%% Integrator
INTEGRATOR = {'CSE_B','GGL_std','GGL_VI'};
DT = [0.05 0.01 0.005 0.001 0.0005 0.0001];
INT_PARA = [NaN NaN; NaN NaN; NaN NaN];
T_0   = 0;
T_END = 1;

%% Solver Method
MAX_ITERATIONS = 40;
TOLERANCE      = 1E-09;

%% Postprocessing
shouldAnimate   = false;
plot_quantities = {};
should_export         = false;
should_export_figures = false;
export_path           = 'scratch/';

%% Write variables into a .mat-File
save(mfilename);