function CONFIG =  metis_load(INPUT_FILE)
%% Function: sets up the METIS workspace and loads all configurations
% Author: Philipp
% date: 02.12.2020

% Define all the parameters that METIS needs to run a simulation in a
% .m-file
run(INPUT_FILE);

% Load input variables into CONFIG-struct and delete unnecessary .mat-File
CONFIG = load([INPUT_FILE,'.mat']);
delete *.mat

end
