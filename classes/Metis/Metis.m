%% Class: Metis
% Conducts all the computations, stores Input-Parameters, sets up
% variables, initialises classes and solving methods
classdef Metis

    properties

        %% Computation parameters
        % Have to be given in an input-file one-to-one. Leaving out one or
        % another will immediately lead to wrong or non-existing results.
        DIM
        DT
        EXT_ACC
        ALL_INTEGRATOR
        INTEGRATOR
        INT_PARA
        INDI_VELO
        MASS
        MAX_ITERATIONS
        Q_0
        SYSTEM
        TOLERANCE
        T_0
        T_END
        V_0
        optime
        NUM_ITER
        MEAN_ITER
        ALPHA_0

        %% Postprocessing parameters
        % Can be given in an input-file (not necessary for computing)
        plot_quantities
        shouldAnimate
        should_export
        export_path
        should_export_figures
        matlab2tikz_directory
        log_file_ID

        %% Solution quantities
        % Will be filled by Metis
        t
        z
        H
        T
        V
        E
        D
        L
        Hdiff
        Ediff
        Ldiff
        diss_work
        constraint_position
        constraint_velocity
        constraint_forces
        external_torque
        alpha_from_q
        mixed_quantity_difference

        %% Error analysis parameters
        % Only necessary if you want to conduct an error analysis for
        % several timestep-sizes
        ALL_DT
        CONV_QUANTITY
        matrix_error_analysis

    end

    methods

        %% Constructor: sets up the METIS workspace and loads all configurations
        function [self, this_system, this_integrator, this_solver] = Metis(INPUT_FILE, num_dt, num_int)

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

        %% Function: Create log-file
        function self = create_log_file(self)

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

        %% Function: Check if user input is valid
        function check_user_input(self)

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

                    %% Function: Check if a class is available
                    function is_avaiable = is_class_available(~, directory, this_class)

                        % Get all present class-names
                        valid_classes = dir(fullfile(directory, '*.m'));

                        % Security check for folders
                        isfile = ~[valid_classes.isdir];
                        filenames = {valid_classes(isfile).name};

                        % Check whether this_class is an available file in the
                        % directory
                        is_avaiable = any(strcmp(filenames, [this_class, '.m']));

                    end

                    %% Function: reads input-file and stores it in the attributes
                    function self = get_config_input(self, INPUT_FILE, num_dt, num_int)

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

                    %% Function: Defines system, integrator and solver objects
                    function [self, this_system, this_integrator, this_solver] = define_classes(self)

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

                    %% Function: set up the matrix with simulation results
                    function self = set_solution_matrix(self, this_integrator, this_system)

                        % Allocate space
                        self.z = zeros(this_integrator.NT+1, this_integrator.nVARS);

                        % Set initial values and time-variable
                        self.z(1, :) = this_integrator.set_initial_condition(self, this_system);
                        self.t = this_integrator.t;

                    end


                end

            end
