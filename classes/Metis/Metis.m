classdef Metis 
    %% Class: Metis
    % Conducts all the computations, stores Input-Parameters, sets up 
    % variables, initilises classes and solving methods    
    
    properties
        %% Computation parameters
        % Have to be given in an input-file one-to-one. Leaving out one or
        % another will immediately lead to wrong or non-existing results.
        DIM
        DT
        EXT_ACC
        INTEGRATOR
        MASS
        MAX_ITERATIONS
        Q_0
        SYSTEM
        TOLERANCE
        T_0
        T_END
        V_0
        optime
        
        %% Postprocessing parameters
        % Can be given in an input-file (not necessary for computing)
        plot_quantities
        shouldAnimate
        should_export
        export_path
        should_export_figures
        
        %% Solution quantities 
        % Will be filled by Metis
        t
        z
        H
        T
        V
        J
        Hdiff
        Jdiff
        constraint_position
        constraint_velocity
        
        %% Error analysis parameters
        % Only necessary if you want to conduct an error analysis for
        % several timestep-sizes
        ALL_DT
        
    end
    
    methods
        
        function [self,this_system, this_integrator, this_solver] = Metis(INPUT_FILE)
            %% Constructor: sets up the METIS workspace and loads all configurations
            % Author: Philipp
            % date: 02.12.2020
            
            % Clear workspace, close all windows and clear command window 
            close all;clc;

            fprintf('**************************************************** \n');
            fprintf('  \n');
            fprintf(' METIS - Computing constrained mechanical systems  \n');
            fprintf(' ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  \n');
            fprintf('     by: Philipp Kinon, B.Sc. \n');
            fprintf('         Institute of Mechanics (IFM) \n')
            fprintf('         Karlsruhe Institute of Technology (KIT) \n')
            fprintf('  \n');
            fprintf('**************************************************** \n');
            fprintf('  \n');
            
            %% Set attributes from config file           
            self = self.get_config_input(INPUT_FILE);
            
            %% Check if user input is valid
            self.check_user_input();
            
            %% Define other classes
            [this_system, this_integrator, this_solver] = self.define_classes();
            
            %% Solution Matrix to store results (inc. IVs)          
            self = self.set_solution_matrix(this_integrator,this_system);
            
        end
        
        function check_user_input(self)
            %% Function: Check if user input is valid
            
            %% Check completeness of input
            % Check if all necessary input quantities are given by user
            necessary_input_list = {'DIM','DT','EXT_ACC','INTEGRATOR','MASS',...
                                    'MAX_ITERATIONS','Q_0','SYSTEM','TOLERANCE',...
                                    'T_0','T_END','V_0'};
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
                fprintf('     Complete user input.                 \n');
                fprintf('  \n');
                fprintf('**************************************************** \n');
                fprintf('  \n');
            end
            
            %% Check validity of given integrator and system        
            % Check if user-input is available in class directories
            is_correct_integrator = self.is_class_available('classes/Integrator',self.INTEGRATOR);
            is_correct_system     = self.is_class_available('classes/System',self.SYSTEM);

            if ~is_correct_integrator 
                error('User input for integrator not available.');
            elseif ~is_correct_system 
                error('User input for system not available.');
            else
                fprintf(['     System: ',self.SYSTEM,' \n']);
                fprintf(['     Integrator: ',self.INTEGRATOR,'\n']);
                fprintf('  \n');
                fprintf('**************************************************** \n');
                
            end

        end
        
        function is_avaiable = is_class_available(~,directory,this_class)
            %% Function: Check if a class is available
            
            % Get all present class-names
            valid_classes     = dir(fullfile(directory,'*.m'));
            
            % Security check for folders
            isfile    = ~[valid_classes.isdir];
            filenames = {valid_classes(isfile).name};
            
            % Check whether this_class is an available file in the
            % directory
            is_avaiable = any(strcmp(filenames,[this_class,'.m']));
            
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
        
        function self = set_solution_matrix(self,this_integrator,this_system)
           
            self.z       = zeros(this_integrator.NT, this_integrator.nVARS);
            z0           = this_integrator.set_initial_condition(self,this_system);
            self.z(1, :) = z0;
            self.t       = this_integrator.t;
            
        end
        
        
    end
    
end