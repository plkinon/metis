classdef Metis 
    %% Class: Metis
    % Conducts all the computations, stores Input-Parameters, sets up 
    % variables, initilises classes and solving methods    
    
    properties
        % All the necessary parameters, initial values and methods used.
        % Have to be given in an input-file.
        DIM
        DT
        EXT_ACC
        INTEGRATOR
        MASS
        MAX_ITERATIONS
        Q_0
        SOLVER
        SYSTEM
        TOLERANCE
        T_0
        t
        T_END
        V_0
        plot_quantities
        shouldAnimate
        ALL_DT
        z
        H
        T
        V
        J
        Hdiff
        Jdiff
        constraint_position
        constraint_velocity
        
    end
    
    methods
        
        function [self,this_system, this_integrator, this_solver] = Metis(INPUT_FILE)
            %% Constructor: sets up the METIS workspace and loads all configurations
            % Author: Philipp
            % date: 02.12.2020
            
            % Clear workspace, close all windows and clear command window 
            close all;clc;

            fprintf('************************************************** \n');
            fprintf(' METIS - Computing constrained mechanical systems \n');
            fprintf('************************************************** \n');
            
            %% Set attributes from config file           
            self = self.get_config_input(INPUT_FILE);
            
            %% Check if user input is valid
            self.check_user_input();
            
            %% Define other classes
            [this_system, this_integrator, this_solver] = self.define_classes();
            
            %% Solution Matrix to store results
            self.z       = zeros(this_integrator.NT, this_integrator.nVARS);
            self.z(1, :) = [self.Q_0', (this_system.MASS_MAT * self.V_0)', this_integrator.LM0'];
            self.t       = this_integrator.t;
            
        end
        
        function check_user_input(self)
            %% Check if user input is valid
            
            % Set strings for directories that contain class m-files
            directory_integrator = 'classes/Integrator';
            directory_system     = 'classes/System';
            
            % Check if user-input is available
            is_correct_integrator = self.is_class_available(directory_integrator,self.INTEGRATOR);
            is_correct_system     = self.is_class_available(directory_system,self.SYSTEM);

            if ~is_correct_integrator 
                error('User input for integrator not available.');
            elseif ~is_correct_system 
                error('User input for system not available.');
            else
                fprintf('               Valid user input.                   \n');
                fprintf('************************************************** \n');
            end

        end
        
        function state = is_class_available(~,directory,this_class)
            %% Function: Check if a class is available
            
            % Get all present class-names
            valid_classes     = dir(fullfile(directory,'*.m'));
            
            % Security check for folders
            isfile    = ~[valid_classes.isdir];
            filenames = {valid_classes(isfile).name};
            
            % Default
            state = false;
            
            % Check whether this_class is an available file in the
            % directory
            if any(strcmp(filenames,[this_class,'.m']))
                state = true;
            end
            
        end
        
        function self = get_config_input(self,INPUT_FILE)
            %% Function: reads input-file and stores it in the attributes
            
            % Define all the parameters that METIS needs to run a simulation in a
            % .m-file
            run(INPUT_FILE);

            % Load input variables into CONFIG-struct and delete unnecessary .mat-File
            configstruct = load([INPUT_FILE,'.mat']);
            delete *.mat
            configstruct = rmfield(configstruct,{'self','INPUT_FILE'});

            % Converts structure s to an object of class classname.
            for fn = fieldnames(configstruct)'    %enumerat fields
               try
                   self.(fn{1}) = configstruct.(fn{1});   %and copy
               catch
                   error('Unknown config parameters given: %s', fn{1});
               end
            end
            
        end
        
        function [this_system, this_integrator, this_solver] = define_classes(self)
            %% Function: Defines problem, integrator and solver objects
            % Author: Philipp Kinon
            % date: 03.12.2020

            %% System
            % Takes user-defined string and evaluates it as the constructor of a class
            % system
            this_system = feval(self.SYSTEM,self);

            %% Integrator
            % Define Integrator from Class (same procedure as for system)
            this_integrator = feval(self.INTEGRATOR,self,this_system);
                     
            %% Solver
            % Define Solver from class
            this_solver = Solver(self);

        end
        
        
    end
    
end