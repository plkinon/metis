%% Class: Solver that conducts the computation
% Time-stepping, Newton Rhapson method and numerical tangent if needed
classdef Solver

    properties
        MAX_ITERATIONS
        TOLERANCE
    end

    methods

        %% Function: Initialise the solver
        function self = Solver(this_simulation)

            % Set Newton tolerance and max. iterations from simulation
            self.TOLERANCE = this_simulation.TOLERANCE;
            self.MAX_ITERATIONS = this_simulation.MAX_ITERATIONS;

        end

        %% Function: Procedure of time-stepping
        function this_simulation = solve(self, this_simulation, this_system, this_integrator)

            % Read solution vector and set present zn1 to initial value
            z = this_simulation.z;
            zn1 = z(1, :)';

            % Command window output
            fprintf('  \n');
            fprintf('     Starting computation ... \n');
            fprintf('  \n');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Starting computation ... ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');

            % analysis of computation time
            n_eval = 1; % number of computations
            stop_time = zeros(n_eval, 1); % allocate space for simulation times

            for j = 1:n_eval

                % Start to count time
                init_time = tic();
                this_simulation.NUM_ITER = zeros(this_integrator.NT,1);
                
                % Actual time-stepping
                for i = 1:this_integrator.NT

                    % Update solution vector from previous iteration
                    zn = zn1;

                    % Initial Guess
                    zn1 = self.newton_initial_guess(zn1,this_system,this_simulation, this_integrator.DT);

                    % Print current time-step
                    fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'----------------------------------------------------');
                    fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),['     Time: t = ',num2str(i*this_integrator.DT),', dt = ', num2str(this_integrator.DT),', step no. ',num2str(i)]);

                    % Conduct iteration
                    [zn1, this_simulation.NUM_ITER(i,1)] = self.newton_iterate(this_integrator, this_system, this_simulation, zn, zn1);

                    % Save update solution for current timestep
                    z(i+1, :) = zn1;

                end

                % Stop counting time and save
                stop_time(j) = toc(init_time);

            end

            % Store solution
            this_simulation.optime = stop_time;
            this_simulation.z = z;
            this_simulation.MEAN_ITER = mean(this_simulation.NUM_ITER);

            % Rearrange unknows if integration schemes requires it
            if ismethod(this_integrator, 'rearrange_unknowns')
                this_simulation.z = this_integrator.rearrange_unknowns(this_simulation, this_system);
            end

            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'---------------------------------------------------- ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Computation suceeded.                           ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
            fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
            fprintf('---------------------------------------------------- \n');
            fprintf('  \n');
            fprintf('     Computation suceeded.                           \n');
            fprintf('  \n');
            fprintf('**************************************************** \n');
            fprintf('  \n');

        end

        %% Funxtion: Initial Guess for Newton-Rhapson-method
        function zn1_guess = newton_initial_guess(~, zn1, this_system, this_simulation, DT)
            
            nDOF = this_system.nDOF;
            M = this_system.MASS_MAT;
            IM = M \ eye(size(M));

            qn1 = zn1(1:nDOF);
            pn1 = zn1(nDOF+1:2*nDOF);

            % Check if integration scheme had independent velocity quantities
            if this_simulation.INDI_VELO == true

                % set velocity vector directly
                vn1 = zn1(2*nDOF+1:3*nDOF);

            else

                % else: compute by means of momenta
                vn1 = IM * pn1;

            end

            qn1_guess = qn1 + DT*vn1;
            zn1_guess = zn1;
            zn1_guess(1:nDOF) = qn1_guess;

        end

        %% Function: Conducts Newton-Rhapson-method to iterate z-vector
        function [zn1, num_iter] = newton_iterate(self, this_integrator, this_system, this_simulation, zn, zn1)

            % Set iteration-index to zero and residual above tolerance
            num_iter = 0;
            residual = self.TOLERANCE * 10;

            % Newton-Rhapson-Method
            while (residual > self.TOLERANCE) && (num_iter <= self.MAX_ITERATIONS)

                % increment iteration index
                num_iter = num_iter + 1;

                % Calculate residual and tangent matrix of present integrator
                % for present system
                [resi, tang] = this_integrator.compute_resi_tang(zn1, zn, this_system);

                % Check if an analytic tangent matrix is implemented
                if isempty(tang)
                    % if not, compute a numerical one
                    tang_num = self.compute_numerical_tangent(this_integrator, this_system, zn1, zn);
                    tang = tang_num;
                end

                % Incrementation of the solution vector
                delta_z = -tang \ resi;
                zn1 = zn1 + delta_z;

                % Compute the residual norm and print current iteration
                residual = max(max(abs(resi)), max(abs(delta_z)));
                fprintf(this_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),['     Iteration ',num2str(num_iter),', residual = ',num2str(residual)]);

            end

        end

        %% Function: numerical tangent
        % Computes numerical tangent matrix for a given residual defined by
        % integrator and system zn1 and zn
        function tang_num = compute_numerical_tangent(~, this_integrator, this_system, zn1, zn)

            % Pre-allocate tangent matrix, initialized with zeros:
            N = length(zn1);
            tang_num = zeros(N, N);

            % Define epsilon which manipulates solution vector
            epsilon = 1e-10;

            % For-loop which computes column by column of numerical tangent matrix
            for j = 1:N

                % Save current entry of solution vector
                zsave = zn1(j);

                % Increment the jth component of the solution vector
                delp = epsilon * (1.0 + abs(zn1(j)));
                zn1(j) = zsave + delp;
                [R1, ~] = this_integrator.compute_resi_tang(zn1, zn, this_system);

                % Decrement the jth component of the solution vector
                zn1(j) = zsave - delp;
                [R2, ~] = this_integrator.compute_resi_tang(zn1, zn, this_system);

                % Restore the original vector
                zn1(j) = zsave;

                % Compute the approximate tangent matrix entry
                tang_num(:, j) = (R1 - R2) / (2 * delp);

            end

        end

    end

end