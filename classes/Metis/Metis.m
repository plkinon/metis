classdef Metis 
 
    properties
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
                   warning('Could not copy field %s', fn{1});
               end
            end
            
        end
        
        function [this_system, this_integrator, this_solver] = initialise(self)
            % Initialises problem, integrator and solver objects
            % Author: Philipp Kinon
            % date: 03.12.2020

            % Clear workspace, close all windows and clear command window 
            close all; clc;

            fprintf('************************************************** \n');
            fprintf(' METIS - Computing constrained mechanical systems \n');
            fprintf('************************************************** \n');

            %% Check if user input is valid
            check_user_input(self);

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