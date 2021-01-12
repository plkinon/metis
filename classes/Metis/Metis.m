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
        T_END
        V_0
        plot_quantities
        shouldAnimate
        INPUT_FILE
        ALL_DT
        
    end
    
    methods
        
        function self = Metis(INPUT_FILE)
            
            % Constructor sets up the METIS workspace and loads all configurations
            % Author: Philipp
            % date: 02.12.2020
            
            % Clear workspace, close all windows and clear command window 
            close all;clc;

            fprintf('************************************************** \n');
            fprintf(' METIS - Computing constrained mechanical systems \n');
            fprintf('************************************************** \n');
            
            % Define all the parameters that METIS needs to run a simulation in a
            % .m-file
            run(INPUT_FILE);

            % Load input variables into CONFIG-struct and delete unnecessary .mat-File
            configstruct = load([INPUT_FILE,'.mat']);
            delete *.mat
            configstruct = rmfield(configstruct,'self');
            
            % Converts structure s to an object of class classname.
            for fn = fieldnames(configstruct)'    %enumerat fields
               try
                   self.(fn{1}) = configstruct.(fn{1});   %and copy
               catch
                   error('Unknown config parameters given: %s', fn{1});
               end
            end
            
        end
        
        function check_user_input(self)
            %% Check if user input is valid
            
            % Set strings for directories that contain class m-files
            directory_integrator = 'classes/Integrator';
            directory_system     = 'classes/System';
            directory_solver     = 'classes/Solver';
            
            % Check if user-input is available
            is_correct_integrator = self.is_class_available(directory_integrator,self.INTEGRATOR);
            is_correct_system     = self.is_class_available(directory_system,self.SYSTEM);
            is_correct_solver     = self.is_class_available(directory_solver,self.SOLVER);

            if ~is_correct_integrator 
                error('User input for integrator not available.');
            elseif ~is_correct_system 
                error('User input for system not available.');
            elseif ~is_correct_solver
                error('User input for solver not available.');
            else
                fprintf('               Valid user input.                   \n');
                fprintf('************************************************** \n');
            end

        end
        
        function state = is_class_available(~,directory,this_class)
            
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
        
        
        function [this_system, this_integrator, this_solver] = initialise(self)
            % Initialises problem, integrator and solver objects
            % Author: Philipp Kinon
            % date: 03.12.2020

            %% Check if user input is valid
            self.check_user_input();

            %% System
            % Takes user-defined string and evaluates it as the constructor of a class
            % system
            this_system = feval(self.SYSTEM,self);

            %% Integrator
            % Define Integrator from Class (same procedure as for system)
            this_integrator = feval(self.INTEGRATOR,self);

            %% Apply integrator to system and vice versa
            % Initialise integrator for current system
            this_integrator = this_integrator.initialise(self,this_system);

            % Intialise System for current integrator
            this_system = this_system.initialise(self,this_integrator);

            %% Solver
            % Define Solver from class
            this_solver = feval(self.SOLVER,self);

        end
        
        
    end
    
end