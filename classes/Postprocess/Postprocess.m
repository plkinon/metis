classdef Postprocess

    methods

        function self = Postprocess()
            % Function to initialise the postprocessing
            set(0, 'defaultfigureposition', [850, 450, 1200, 600])
        end

        function animation(~, this_problem, this_integrator, CONFIG)
            % Function which displays an animation of the trajectory if
            % desired by user

            if CONFIG.shouldAnimate == true

                fig = figure();
                DIM = this_problem.DIM;
                q = this_problem.z(:, 1:DIM);
                l = this_problem.GEOM(1);
                NT = this_integrator.NT;

                title(strcat(this_integrator.NAME, ': Trajectory'))
                axis equal
                axis([-1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l]);
                xlabel('x');
                ylabel('y');
                zlabel('z');
                grid on;

                xa = q(1, 1);
                ya = q(1, 2);
                if DIM == 3
                    za = q(1, 3);
                else
                    za = 0;
                end

                for j = 1:NT + 1

                    cla(fig);
                    hold on

                    %% Current position
                    x = q(j, 1);
                    y = q(j, 2);
                    if DIM == 3
                        z = q(j, 3);
                    else
                        z = 0;
                    end

                    %% Reference sphere
                    plot3(xa, ya, za, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    hold on

                    %% Reference constraint
                    x3 = [0; xa];
                    y3 = [0; ya];
                    z3 = [0; za];
                    plot3(x3, y3, z3, 'k', 'LineWidth', 1);

                    %% current position of the mass
                    hold on
                    if DIM == 3
                        plot3(q(1:j, 1), q(1:j, 2), q(1:j, 3), 'k');
                    else
                        plot3(q(1:j, 1), q(1:j, 2), zeros(j, 1), 'k');
                    end
                    plot3(x, y, z, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    grid on

                    %% current position of the constraint
                    x3 = [0; x];
                    y3 = [0; y];
                    z3 = [0; z];
                    plot3(x3, y3, z3, 'k', 'linewidth', 1);
                    if DIM == 2
                        view(0, 90)
                    else
                        view(136, 23)
                    end

                    drawnow
                    
                end

            end

        end

        function this_problem = compute(~, this_problem)
            % Function which computes energetic quantities, ang. Mom. , ...
            % as functions of time for a given q and p-vector

            DIM = this_problem.DIM;
            q = this_problem.z(:, 1:DIM);
            p = this_problem.z(:, DIM+1:2*DIM);
            NT = size(q, 1);
            M = this_problem.MASS_MAT;
            IM = M\eye(size(M));

            T = zeros(NT, 1);
            V = zeros(NT, 1);
            H = zeros(NT, 1);
            J = zeros(NT, 3);
            diffH = zeros(NT-1, 1);
            diffJ = zeros(NT-1, 3);
            constraint_position = zeros(NT,1);
            constraint_velocity = zeros(NT,1);

            for j = 1:NT
                T(j) = 1 / 2 * p(j, :) * M * p(j, :)';
                V(j) = this_problem.potential(q(j,:)');
                H(j) = T(j) + V(j);
                constraint_position(j) = this_problem.constraint(q(j,:)');
                constraint_velocity(j) = this_problem.constraint_gradient(q(j,:)')*IM*p(j,:)';

                if DIM == 3
                    J(j, :) = cross(q(j, :), p(j, :));
                end

                if j > 1
                    diffH(j-1) = H(j) - H(j-1);
                    diffJ(j-1, :) = J(j, :) - J(j-1, :);
                end

            end

            this_problem.T = T;
            this_problem.V = V;
            this_problem.H = H;
            this_problem.J = J;
            this_problem.Hdiff = diffH;
            this_problem.Jdiff = diffJ;
            this_problem.constraint_position = constraint_position;
            this_problem.constraint_velocity = constraint_velocity;

        end

        function plot(~,CONFIG, this_problem, this_integrator)
                      
            H = this_problem.H;
            V = this_problem.V;
            T = this_problem.T;
            t = this_integrator.t;
            J = this_problem.J;
            diffH = this_problem.Hdiff;
            diffJ = this_problem.Jdiff;
            g_pos = this_problem.constraint_position;
            g_vel = this_problem.constraint_velocity;
            
            Mmin = min([min(V), min(T), min(H)]);
            Mmax = max([max(V), max(T), max(H)]);
            
            for i = 1:length(CONFIG.plot_quantities)
                
                figure();
                grid on;
                quantity = CONFIG.plot_quantities{i};
            
                switch quantity

                    case 'energy'

                        %plots Energy quantities over time
                        plotline = plot(t, T, t, V, t, H);
                        ylim([Mmin - 0.1 * abs(Mmin), 1.1 * Mmax])
                        title(strcat(this_integrator.NAME, ': Energy'))
                        legend('T', 'V', 'H')
                        xlabel('t');
                        ylabel('Energy');

                    case 'angular_momentum_3'

                        %plots the ang. Mom. about dim-axis over time
                        axis = 3;
                        plotline = plot(t, J(:,axis));
                        title(strcat(this_integrator.NAME, ': Angular Momentum'))
                        xlabel('t');
                        ylabel(strcat('Angular Momentum about {} ', num2str(axis), '-axis'));

                    case 'energy_difference'

                        %plots the increments of Hamiltonian over time
                        plotline = plot(t(2:end), diffH);
                        title(strcat(this_integrator.NAME, ': Energy difference'))
                        xlabel('t');
                        ylabel('H_{n+1}-H_{n}');
                        legend(strcat('std(H)=', num2str(std(H))));

                    case 'angular_momentum_difference_3'

                        %plots the increments of ang. Mom. over time about desired axis (dim)
                        axis = 3;             
                        plotline = plot(t(2:end), diffJ(:,axis));
                        title(strcat(this_integrator.NAME, ': Ang. mom. - difference'));
                        xlabel('t');
                        ylabel(strcat('J_{n+1}-J_{n} about {} ', num2str(axis), '-axis'));
                        legend(strcat('std(J)=', num2str(std(J(:,axis)))));


                    case 'constraint_position'

                        %plots the position constraint and their violations by integration
                        plotline = plot(t,g_pos);
                        title(strcat(this_integrator.NAME, ': Constraint on position level'));
                        xlabel('t');
                        ylabel('g(t)');

                    case 'constraint_velocity'

                        %plots the velocity constraint and their violations by integration
                        plotline = plot(t,g_vel);
                        title(strcat(this_integrator.NAME, ': Constraint on velocity level'));
                        xlabel('t');
                        ylabel('G*v(t)');    

                    otherwise

                        error(strcat('No plotting routine for ',quantity,'  this quantity defined'));

                end
            
                set(plotline, 'linewidth', 1.5);
                grid on;
                
            end
        
        end
        
    end

end
