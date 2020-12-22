classdef Newton

    %% Solver class of Newton Iteration

    properties
        MAX_ITERATIONS
        TOLERANCE
    end

    methods

        function self = Newton(CONFIG)

            %% Initialising the Solver
            self.TOLERANCE  = CONFIG.TOLERANCE;
            self.MAX_ITERATIONS = CONFIG.MAX_ITERATIONS;
            
        end

        function this_problem = solve(self, this_problem, this_integrator)

            %% Procedure of Newton-Rhapson-Method for time-stepping
            z = this_problem.z;
            zn1 = z(1, :)';

            for i = 1:this_integrator.NT

                %% Time-Stepping

                zn = zn1;
                k = 0;
                residual = 1e05;

                fprintf('------------------------------------------\n');
                fprintf('Time: t = %.3f, dt = %.3f, step no. %.0f \n', i*this_integrator.DT, this_integrator.DT, i);

                while residual > self.TOLERANCE && k <= self.MAX_ITERATIONS

                    %% Newton-Rhapson-Method
                    k = k + 1;

                    % Calculate residual and tangent matrix
                    [resi,tang] = this_integrator.compute_resi_tang(zn1, zn, this_problem);
                    
                    % Check if an analytic tangent matrix is implemented,
                    % if not, compute a numerical one
                    if isempty(tang)
                        tang = self.compute_numerical_tangent(this_integrator, this_problem, zn1, zn);
                    end
                                          
                    %% Incrementation of the solution vector
                    dz  = -tang \ resi;
                    zn1 = zn1 + dz;    

                    residual = max(max(abs(resi)), max(dz));

                    fprintf('Iteration %.0f) residual = %.4d \n', k, residual);


                end

                z(i+1, :) = zn1;

            end

            this_problem.z = z;
        end
        
        function tang_num = compute_numerical_tangent(~, this_integrator,this_problem,zn1, zn)
            
            %% FUNCTION: numerical tangent
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