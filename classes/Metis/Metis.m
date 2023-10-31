%% Class: Metis
% Conducts all the computations, stores Input-Parameters, sets up
% variables, initialises classes and solving methods
classdef Metis
% Conducts all the computations, stores Input-Parameters, sets up
% variables, initialises classes and solving methods:
% Computation parameters have to be given in an input-file one-to-one. Leaving out one or
% another will immediately lead to wrong or non-existing results.
% Postprocessing parameters can be given in an input-file (not necessary for computing).
% Solution quantities will be filled by Metis.
% Error analysis parameters are only necessary if you want to conduct an error analysis for 
% several timestep-sizes.

    properties

        DIM % number of spatial dimensions
        DT % fixed time step size
        EXT_ACC % external acceleration vector
        INTEGRATOR % name of desired integrator
        INT_PARA % parameters of integrator  
        INDI_VELO % Boolean if the integrators features independent velocities
        MASS % masses of the system
        MAX_ITERATIONS % maximal number of iterations for the Newton's method
        Q_0 % initial configuration
        SYSTEM % name of the desired system
        TOLERANCE % tolerance of Newton's method
        T_0 % start time
        T_END % end time
        V_0 % initial velocity
        optime % saves time which was needed for the computation (in seconds)
        NUM_ITER % saves number of iterations
        MEAN_ITER % saves mean value of iterations per time step
        ALPHA_0 % initial value for mixed variables (if present)

    
        plot_quantities % determines which quantities will be plotted in the postprocessing
        shouldAnimate % determines if there is an animation
        should_export % determines if data should be exported
        export_path % data path for export
        should_export_figures % determines if figures should be exported
        matlab2tikz_directory % path for the matlab2tikz extension
        log_file_ID % ID of log file

 
        t %time
        z %solution vector
        H % total energy
        T %kinetic energy
        V %potential energy
        E %generalized energy function
        D %dissipation work
        L % angular momentum
        Hdiff % incremental difference of total energy
        Ediff % incremental difference of generalized energy function
        Ldiff % incremental difference of angular momentum
        diss_work % accumulated dissipated wprk
        constraint_position % value of constraint on position level
        constraint_velocity%  value of constraint on velocity/momentum level
        constraint_forces % constraint forces
        external_torque % external torques
        alpha_from_q % determine mixed quantity from configuration
        mixed_quantity_difference % % incremental difference of mixed quantity

        
        ALL_DT % all time step sizes for a convergence analysis
        ALL_INTEGRATOR % all integrators for a convergence analysis
        CONV_QUANTITY % determines which quantity should be analyzed in the convergence analysis

    end

    methods

        function [self, this_system, this_integrator, this_solver] = Metis(INPUT_FILE, num_dt, num_int)
        % Constructor: sets up the METIS workspace and loads all configurations
            
            % Clear workspace, close all windows and clear command window
            close all;
            clc;

            fprintf('**************************************************** \n');
            fprintf('  \n');
            fprintf(' METIS - Computing constrained mechanical systems  \n');
            fprintf(' ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  \n');
            fprintf('     by: Philipp Kinon, M.Sc. \n');
            fprintf('         Institute of Mechanics (IFM) \n');
            fprintf('         Karlsruhe Institute of Technology (KIT) \n');
            fprintf('  \n');
            fprintf('**************************************************** \n');
            fprintf('  \n');

            %Create log-file
            self = self.create_log_file();

            % Start visual output to command window
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),' METIS - Computing constrained mechanical systems  ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),' ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'     by: Philipp Kinon, M.Sc. ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'         Institute of Mechanics (IFM) ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'         Karlsruhe Institute of Technology (KIT) ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),' ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
            fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');

            % Set attributes from config file
            self = self.get_config_input(INPUT_FILE, num_dt, num_int);

            % Check if user input is valid
            self.check_user_input();

            % Define other classes
            [self, this_system, this_integrator, this_solver] = self.define_classes();

            % Solution Matrix to store results (inc. IVs)
            self = self.set_solution_matrix(this_integrator, this_system);

        end


        function self = create_log_file(self)
        % Create log-file

            % Create new log file
                if ~exist('metis.log', 'file')

                    self.log_file_ID = fopen(fullfile(self.export_path, 'metis.log'), 'a');

                else

                    %if there already is one, overwrite it
                    fprintf('     Overwriting existing log file.            \n');
                    fprintf('  \n');
                    fprintf('**************************************************** \n');
                    fprintf('  \n');
                    warning('off')
                    delete('metis.log');
                    self.log_file_ID = fopen(fullfile(self.export_path, 'metis.log'), 'a');
                    warning('on')

                end

            % Start logging into the log file
            self.log_file_ID = fopen(fullfile(self.export_path, 'metis.log'), 'a');
            % Check if log-file can be opened/edited
            if self.log_file_ID == -1
              error('Cannot open log file.');
            end

        end

        function check_user_input(self)
        % Check if user input is valid

            % Check if all necessary input quantities are given by user
            necessary_input_list = {'DIM', 'DT', 'EXT_ACC', 'INTEGRATOR', 'MASS', 'MAX_ITERATIONS', 'Q_0', 'SYSTEM', 'TOLERANCE', 'T_0', 'T_END', 'V_0'};
            is_user_input_complete = true;

            for i = 1:length(necessary_input_list)
                if isempty(self.(necessary_input_list{i}))
                    is_user_input_complete = false;
                end
            end

            % Throw error if necessary input quantity is missing
            if ~is_user_input_complete
                error('User input is incomplete. Check for missing input values.');
            else
                fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'     Complete user input.                           ');
                fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'****************************************************');
                fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf('     Complete user input.                            \n');
                fprintf('  \n');
                fprintf('**************************************************** \n');
                fprintf('  \n');
            end

                % Check validity of given integrator and system
                is_correct_integrator = self.is_class_available('classes/Integrator', self.INTEGRATOR);
                is_correct_system = self.is_class_available('classes/System', self.SYSTEM);

                if ~is_correct_integrator
                    error('User input for integrator not available.');
                elseif ~is_correct_system
                    error('User input for system not available.');
                else
                    fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),['     System: ', self.SYSTEM, ' ']);
                    fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),['     Integrator: ', self.INTEGRATOR, '']);
                    fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                    fprintf(self.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
                    fprintf(['     System: ', self.SYSTEM, ' \n']);
                    fprintf(['     Integrator: ', self.INTEGRATOR, '\n']);
                    fprintf('  \n');
                    fprintf('**************************************************** \n');

                end

        end

                    function is_avaiable = is_class_available(~, directory, this_class)
                    % Check if a class is available

                        % Get all present class-names
                        valid_classes = dir(fullfile(directory, '*.m'));

                        % Security check for folders
                        isfile = ~[valid_classes.isdir];
                        filenames = {valid_classes(isfile).name};

                        % Check whether this_class is an available file in the
                        % directory
                        is_avaiable = any(strcmp(filenames, [this_class, '.m']));

                    end

                    function self = get_config_input(self, INPUT_FILE, num_dt, num_int)
                    % reads input-file and stores it in the attributes
                        % Define all the parameters that METIS needs to run a
                        % simulation in a .m-file
                        run(INPUT_FILE);

                        % Load input variables into CONFIG-struct and delete unnecessary .mat-File
                        configstruct = load([INPUT_FILE, '.mat']);
                        delete([INPUT_FILE, '.mat']);
                        configstruct = rmfield(configstruct, {'self', 'INPUT_FILE', 'num_dt', 'num_int'});

                        % Converts structure s to an object of class classname.
                        for fn = fieldnames(configstruct)' %enumerat fields
                            try
                                self.(fn{1}) = configstruct.(fn{1}); %and copy
                            catch
                                error('Unknown config parameters given: %s', fn{1});
                            end
                        end

                        % For multiple simulations at once set current time step size
                        self.ALL_DT = self.DT;
                        self.DT = self.DT(num_dt);

                        % For multiple simulations at once set current integrator
                        self.ALL_INTEGRATOR = self.INTEGRATOR;
                        if ~ischar(self.INTEGRATOR)
                            self.INTEGRATOR = self.INTEGRATOR{num_int};
                        end

                        % For multiple simulations at once set current integrator parameters
                        self.INT_PARA = self.INT_PARA(num_int, :);

                    end

                    function [self, this_system, this_integrator, this_solver] = define_classes(self)
                    % Defines system, integrator and solver objects
                        % System: Takes user-defined string and evaluates it as the
                        % constructor of a class system
                        this_system = feval(self.SYSTEM, self);

                        % Integrator: Define from Class (same procedure as for system)
                        this_integrator = feval(self.INTEGRATOR, self, this_system);
                        % Clear given integrator parameter if integrator has none
                        if ~this_integrator.hasPARA
                            clear this_simulation.INT_PARA
                        end
                        self.INDI_VELO = this_integrator.INDI_VELO;

                        % Solver: Define from class
                        this_solver = Solver(self);

                    end

                    
                    function self = set_solution_matrix(self, this_integrator, this_system)
                    %% set up the matrix with simulation results

                        % Allocate space
                        self.z = zeros(this_integrator.NT+1, this_integrator.nVARS);

                        % Set initial values and time-variable
                        self.z(1, :) = this_integrator.set_initial_condition(self, this_system);
                        self.t = this_integrator.t;

                    end


                end

            end
