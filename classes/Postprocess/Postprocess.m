classdef Postprocess

    methods

        function self = Postprocess()
            % Function to initialise the postprocessing
            set(0, 'defaultfigureposition', [850, 450, 1200, 600])
        end

        function animation(~, this_problem, this_simulation)
            % Function which displays an animation of the trajectory if
            % desired by user
            
            if this_simulation.shouldAnimate == true
                
                fig = figure();
                this_problem.give_animation(fig,this_simulation);
                title(strcat(this_simulation.INTEGRATOR, ': Trajectory'))

            end

        end

        function this_simulation = compute(~, this_problem,this_simulation)
            % Function which computes energetic quantities, ang. Mom. , ...
            % as functions of time for a given q and p-vector
            nDOF = this_problem.nDOF;
            m = this_problem.mCONSTRAINTS;
            d = this_problem.nBODIES;
            DIM = this_problem.DIM;
            q = this_simulation.z(:, 1:nDOF);
            p = this_simulation.z(:, nDOF+1:2*nDOF);
            NT = size(q, 1);
            M = this_problem.MASS_MAT;
            IM = M\eye(size(M));

            T = zeros(NT, 1);
            V = zeros(NT, 1);
            H = zeros(NT, 1);
            J = zeros(NT, 3);
            diffH = zeros(NT-1, 1);
            diffJ = zeros(NT-1, 3);
            constraint_position = zeros(NT,m);
            constraint_velocity = zeros(NT,m);

            for j = 1:NT
                T(j) = 1 / 2 * p(j, :) * IM' * p(j, :)';
                V(j) = this_problem.internal_potential(q(j,:)') + this_problem.external_potential(q(j,:)');
                H(j) = T(j) + V(j);
                constraint_position(j,:) = this_problem.constraint(q(j,:)')';
                constraint_velocity(j,:) = (this_problem.constraint_gradient(q(j,:)')*IM*p(j,:)')';

                if DIM == 3
                    for k = 1:d
                        J(j, :) = J(j,:) + cross(q(j,(k-1)*DIM+1:k*DIM), p(j, (k-1)*DIM+1:k*DIM));
                    end
                end

                if j > 1
                    diffH(j-1) = H(j) - H(j-1);
                    diffJ(j-1, :) = J(j, :) - J(j-1, :);
                end

            end
                       
            this_simulation.T = T;
            this_simulation.V = V;
            this_simulation.H = H;
            this_simulation.J = J;
            this_simulation.Hdiff = diffH;
            this_simulation.Jdiff = diffJ;
            this_simulation.constraint_position = constraint_position;
            this_simulation.constraint_velocity = constraint_velocity;

        end
        
        function save(~,export_simulation)
            
            % convert timestep-size to string
            DT_string         = strrep(num2str(export_simulation.DT),'.','');
            integrator_string = strrep(export_simulation.INTEGRATOR,'_','');
            % Set export-string-name
            export_folder = [export_simulation.export_path,export_simulation.SYSTEM,'_',integrator_string,'_DT',DT_string,'/'];
                
            
            if export_simulation.should_export
                if ~exist(export_folder, 'dir')
                    mkdir(export_folder)
                else
                    warning('off')
                    rmdir(export_folder,'s')
                    mkdir(export_folder)
                    warning('on')
                    warning('Removed existing output folder.');
                end
                
                % Export
                save([export_folder,'results'],'export_simulation');
                
                if export_simulation.should_export_figures
                    figHandles = findall(0,'Type','figure'); 
                    
                    for i = 1:numel(figHandles)
                        print(figHandles(i), [export_folder,num2str(i)], '-depsc')
                    end
                    
                end
            end
            
        end

        function plot(~,this_simulation)
            
            H = this_simulation.H;
            V = this_simulation.V;
            T = this_simulation.T;
            t = this_simulation.t;
            J = this_simulation.J;
            diffH = this_simulation.Hdiff;
            diffJ = this_simulation.Jdiff;
            g_pos = this_simulation.constraint_position;
            g_vel = this_simulation.constraint_velocity;
            
            Mmin = min([min(V), min(T), min(H)]);
            Mmax = max([max(V), max(T), max(H)]);
            
            for i = 1:length(this_simulation.plot_quantities)
                
                figure();
                grid on;
                quantity = this_simulation.plot_quantities{i};
            
                switch quantity

                    case 'energy'

                        %plots Energy quantities over time
                        plotline = plot(t, T, t, V, t, H);
                        ylim([Mmin - 0.1 * abs(Mmin), 1.1 * Mmax])
                        title(strcat(this_simulation.INTEGRATOR, ': Energy'))
                        legend('T', 'V', 'H')
                        xlabel('t');
                        ylabel('Energy');

                    case 'angular_momentum'

                        %plots the ang. Mom. about dim-axis over time
                        plotline = plot(t, J(:,:));
                        title(strcat(this_simulation.INTEGRATOR, ': Angular Momentum'))
                        legend('J1','J2','J3');
                        xlabel('t');
                        ylabel('Angular Momentum');

                    case 'energy_difference'

                        %plots the increments of Hamiltonian over time
                        plotline = plot(t(2:end), diffH);
                        title(strcat(this_simulation.INTEGRATOR, ': Energy difference'))
                        xlabel('t');
                        ylabel('H_{n+1}-H_{n}');
                        legend(strcat('std(H)=', num2str(std(H))));

                    case 'angular_momentum_difference'

                        %plots the increments of ang. Mom. over time about desired axis (dim)      
                        plotline = plot(t(2:end), diffJ(:,:));
                        title(strcat(this_simulation.INTEGRATOR, ': Ang. mom. - difference'));
                        xlabel('t');
                        legend('J1','J2','J3');
                        ylabel('J_{n+1}-J_{n}');
                        legend(strcat('std(J)=[', num2str(std(J(:,:))),']'));


                    case 'constraint_position'

                        %plots the position constraint and their violations by integration
                        plotline = plot(t,g_pos);
                        title(strcat(this_simulation.INTEGRATOR, ': Constraint on position level'));
                        xlabel('t');
                        ylabel('g(t)');

                    case 'constraint_velocity'

                        %plots the velocity constraint and their violations by integration
                        plotline = plot(t,g_vel);
                        title(strcat(this_simulation.INTEGRATOR, ': Constraint on velocity level'));
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
