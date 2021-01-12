classdef Solver

    %% Solver class that conducts the computation

    properties
        MAX_ITERATIONS
        TOLERANCE
    end

    methods

        function self = Solver(this_simulation)
            %% Initialising the Solver
            self.TOLERANCE  = this_simulation.TOLERANCE;
            self.MAX_ITERATIONS = this_simulation.MAX_ITERATIONS;
            
        end

        function this_problem = solve(self, this_problem, this_integrator)
            %% Procedure of time-stepping
            
            % Read solution vector and set present zn1 to initial value
            z = this_problem.z;
            zn1 = z(1, :)';
            
            %% Time-Stepping
            % Forward time-stepping until end
            for i = 1:this_integrator.NT
      
                % Update solution vector from previous iteration
                zn = zn1;
                                
                % Print current time-step
                fprintf('------------------------------------------\n');
                fprintf('Time: t = %.3f, dt = %.3f, step no. %.0f \n', i*this_integrator.DT, this_integrator.DT, i);
                
                % Conduct iteration
                zn1 = self.newton_iterate(this_integrator, this_problem,zn, zn1);
                
                % Save update solution for current timestep
                z(i+1, :) = zn1;

            end
            
            % Store solution
            this_problem.z = z;
            
        end
        
        function zn1 = newton_iterate(self,this_integrator, this_problem, zn, zn1)
            %% Function: Conducts Newton-Rhapson-method to iterate z-vector
            
            % Set iteration-index to zero and residual above tolerance
            k = 0;
            residual = 1e05;
            
            while residual > self.TOLERANCE && k <= self.MAX_ITERATIONS

                %% Newton-Rhapson-Method
                k = k + 1;
                
                % Calculate residual and tangent matrix of present integrator
                % for present system
                [resi,tang] = this_integrator.compute_resi_tang(zn1, zn, this_problem);

                % Check if an analytic tangent matrix is implemented,
                % if not, compute a numerical one
                if isempty(tang)
                    tang = self.compute_numerical_tangent(this_integrator, this_problem, zn1, zn);
                end

                % Incrementation of the solution vector
                delta_z  = -tang \ resi;
                zn1      = zn1 + delta_z;    

                % Compute the residual norm and print current iteration
                residual = max(max(abs(resi)), max(delta_z));
                fprintf('Iteration %.0f) residual = %.4d \n', k, residual);
                
            end
            
        end
        
        function tang_num = compute_numerical_tangent(~, this_integrator,this_problem,zn1, zn)
            %% Function: numerical tangent
            % Computes numerical tangent matrix for a given residual defined by
            % integrator and problem zn1 and zn

            % Pre-allocate tangent matrix, initialized with zeros:
            N          = length(zn1);
            tang_num   = zeros(N,N);

            % Define epsilon which manipulates solution vector
            epsilon = 1e-10;

            % For-loop which computes column by column of numerical tangent matrix
            for j=1:N

                % Save current entry of solution vector
                zsave=zn1(j);

                % Increment the jth component of the solution vector 
                delp=epsilon*(1.0+abs(zn1(j)));
                zn1(j)=zsave+delp;
                [R1,~] = this_integrator.compute_resi_tang(zn1, zn, this_problem);

                % Decrement the jth component of the solution vector
                zn1(j)=zsave-delp;
                [R2,~] = this_integrator.compute_resi_tang(zn1, zn, this_problem);

                % Restore the original vector 
                zn1(j)=zsave;

                % Compute the approximate tangent matrix entry
                tang_num(:,j)=(R1-R2)/(2*delp);

            end

        end
        
    end
end