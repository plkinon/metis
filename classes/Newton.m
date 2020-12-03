classdef Newton

    %% Solver class of Newton Iteration

    properties
        MAX_ITERATIONS
        TOLERANCE
        NUM_TANGENT
    end

    methods

        function self = Newton(CONFIG)

            %% Initialising the Solver
            self.TOLERANCE  = CONFIG.TOLERANCE;
            self.MAX_ITERATIONS = CONFIG.MAX_ITERATIONS;
            self.NUM_TANGENT = CONFIG.NUM_TANGENT;
            
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
                    if self.NUM_TANGENT == true
                        
                        [resi,~] = this_integrator.compute_resi_tang(zn1, zn, this_problem);
                        tang_num = numerical_tangent(this_integrator, this_problem, zn1, zn);
                        tang     = tang_num;
                        
                    elseif self.NUM_TANGENT == false
                        
                        [resi, tang] = this_integrator.compute_resi_tang(zn1, zn, this_problem);
                        
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
    end
end