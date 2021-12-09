%% Class: Postprocessing
% Distributes animation, computation of postprocessing quantities and plots
classdef Postprocess

    properties

        % Color scheme also for colorblind readers
        color_scheme = {'#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'};

    end


    methods

        %% Function: initialise the postprocessing
        function self = Postprocess()
            % set default position of figure windows
            set(0, 'defaultfigureposition', [850, 450, 1000, 1000])

        end

        %% Function: displays an animation of the trajectory
        function animation(~, this_system, this_simulation)

            if this_simulation.shouldAnimate == true

                fig = figure();
                % accesses routine from system-object
                this_system.give_animation(fig, this_simulation);
                title(strcat(this_simulation.INTEGRATOR, ': Trajectory'))
                current_fig = gcf;
                current_fig.Name = 'final_conf';

            end

        end

        %% Function: computes postprocessing quantities as function of time
        function this_simulation = compute(~, this_system, this_simulation)

            nDOF = this_system.nDOF;
            m = this_system.mCONSTRAINTS;
            DIM = this_system.DIM;
            d = nDOF / DIM;

            NT = int32(this_simulation.T_END/this_simulation.DT) + 1;
            M = this_system.MASS_MAT;
            IM = M \ eye(size(M));

            q = this_simulation.z(:, 1:nDOF);
            p = this_simulation.z(:, nDOF+1:2*nDOF);
            v = zeros(NT, nDOF);

            % Check if integration scheme had independent velocity quantities
            if this_simulation.INDI_VELO == true

                % set velocity vector directly
                v = this_simulation.z(:, 2*nDOF+1:3*nDOF);

            else

                % else: compute by means of momenta
                for j = 1:NT
                    v(j, :) = (IM * p(j, :)')';
                end

            end

            % Allocate space for postprocessing quantities
            T = zeros(NT, 1);
            V = zeros(NT, 1);
            H = zeros(NT, 1);
            L = zeros(NT, 3);
            diffH = zeros(NT-1, 1);
            diffL = zeros(NT-1, 3);
            constraint_position = zeros(NT, m);
            constraint_velocity = zeros(NT, m);

            % Compute quantities for every point in time
            for j = 1:NT

                % Kinetic and potential energy, Hamiltonian
                T(j) = 1 / 2 * v(j, :) * M * v(j, :)';
                %T(j) = 1/2*p(j,:)*IM*p(j,:)';
                V(j) = this_system.internal_potential(q(j, :)') + this_system.external_potential(q(j, :)');
                H(j) = T(j) + V(j);

                % Constraints on position and velocity level
                constraint_position(j, :) = this_system.constraint(q(j, :)')';
                constraint_velocity(j, :) = (this_system.constraint_gradient(q(j, :)') * v(j, :)')';

                if strcmp(this_simulation.INTEGRATOR, 'GGL_VI_mod')
                    % this specific integration scheme relies upon
                    % secondary constraints with momentum formulation and
                    % has non-aligned velocites
                    constraint_velocity(j, :) = (this_system.constraint_gradient(q(j, :)') * IM * p(j, :)')';
                end

                % Compute angular momentum
                if DIM == 3

                    for k = 1:d

                        L(j, :) = L(j, :) + cross(q(j, (k - 1)*DIM+1:k*DIM), p(j, (k - 1)*DIM+1:k*DIM));

                    end

                end

            end

            % Compute time-increments in Hamiltonian and angular momentum
            for j = 1:(NT - 1)

                diffH(j) = H(j+1) - H(j);
                diffL(j, :) = L(j+1, :) - L(j, :);

            end

            % Save computed quantities to simulation-object
            this_simulation.T = T;
            this_simulation.V = V;
            this_simulation.H = H;
            this_simulation.L = L;
            this_simulation.Hdiff = diffH;
            this_simulation.Ldiff = diffL;
            this_simulation.constraint_position = constraint_position;
            this_simulation.constraint_velocity = constraint_velocity;

        end

        %% Function: save results
        function save(~, export_simulation)

            % convert timestep-size to string
            DT_string = strrep(num2str(export_simulation.DT), '.', '');
            integrator_string = strrep(export_simulation.INTEGRATOR, '_', '');

            % Set export-string-name
            export_folder = [export_simulation.export_path, export_simulation.SYSTEM, '_', integrator_string, '_DT', DT_string, '/'];

            % Check if export is desired by user
            if export_simulation.should_export

                % Create new output folder
                if ~exist(export_folder, 'dir')

                    mkdir(export_folder)

                else

                    warning('off')
                    rmdir(export_folder, 's')
                    mkdir(export_folder)
                    warning('on')
                    fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Overwriting existing output folder.           ');
                    fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                    fprintf('     Overwriting existing output folder.            \n');
                    fprintf('  \n');

                end

                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Exporting results.                 ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf('     Exporting results.                 \n');
                fprintf('  \n');

                % Export as .mat file
                save([export_folder, 'results'], 'export_simulation');

                % Check if plots should be exported as well
                if export_simulation.should_export_figures
                    figHandles = findall(0, 'Type', 'figure');

                    % go through existing figures
                    for i = 1:numel(figHandles)

                        % Check if current figure has a name
                        if ~isempty(figHandles(i).Name)
                            export_name = figHandles(i).Name;
                        else
                            export_name = num2str(i);
                        end

                        % Clear current figures title for export
                        set(0, 'currentfigure', figHandles(i));
                        titlestring = get(gca, 'title').String;
                        title('');
                        % Reduce linewidth
                        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);

                        % Export to .eps
                        print(figHandles(i), [export_folder, export_name], '-depsc')
                                        
                        % Export to .tikz
                        warning('off')
                        try
                            % add chosen matlab2tikz directory
                            addpath([export_simulation.matlab2tikz_directory,'/src']);
                            % tries to export via matlab2tikz if available
                            matlab2tikz('figurehandle', figHandles(i), 'height', '\figH', 'width', '\figW', 'filename', [export_folder, export_name, '.tikz'], 'showInfo', false, 'floatformat', '%.7g');
                        catch
                            % if chosen directory is wrong or not existing
                            fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),['     Matlab2Tikz not found at ',export_simulation.matlab2tikz_directory,' .           ']);
                            fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                            fprintf(['     Matlab2Tikz not found at ',export_simulation.matlab2tikz_directory,' .           \n']);
                            fprintf('  \n');
                            % clone current matlab2tikz repository 
                            [~, ~] = system(['git clone git@github.com:matlab2tikz/matlab2tikz.git ',export_simulation.matlab2tikz_directory]);
                            fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),['     Cloning Matlab2Tikz to ',export_simulation.matlab2tikz_directory,'                 ']);
                            fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                            fprintf(['     Cloning Matlab2Tikz to ',export_simulation.matlab2tikz_directory,'                 \n']);
                            fprintf('  \n');
                            addpath([export_simulation.matlab2tikz_directory,'/src']);
                            % Export to .tikz now
                            matlab2tikz('figurehandle', figHandles(i), 'height', '\figH', 'width', '\figW', 'filename', [export_folder, export_name, '.tikz'], 'showInfo', false, 'floatformat', '%.7g');
                        end
                        warning('on')

                        %Reset Title and enlarge linewidth
                        set(0, 'currentfigure', figHandles(i));
                        title(titlestring);
                        set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.5);

                    end

                end

                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Exporting succesful.                 ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
                fprintf('     Exporting succesful.                 \n');
                fprintf('  \n');
                fprintf('**************************************************** \n');
                % Close log-file
                fclose(export_simulation.log_file_ID);

            end

        end

        %% Function: plot specific quantities
        function plot(self, this_simulation)
            % Function for plotting postprocessing results

            H = this_simulation.H;
            V = this_simulation.V;
            T = this_simulation.T;
            t = this_simulation.t;
            L = this_simulation.L;
            diffH = this_simulation.Hdiff;
            diffL = this_simulation.Ldiff;
            g_pos = this_simulation.constraint_position;
            g_vel = this_simulation.constraint_velocity;

            % go through desired output quantities
            for i = 1:length(this_simulation.plot_quantities)

                % get current quantity
                fig = figure();
                grid on;
                quantity = this_simulation.plot_quantities{i};
                integrator_string = strrep(this_simulation.INTEGRATOR, '_', '-');

                switch quantity

                    case 'energy'

                        %plots Energy quantities over time
                        if isnan(H(end))
                            plotline = plot(t(1:end-1), T(1:end-1), t(1:end-1), V(1:end-1), t(1:end-1), H(1:end-1));
                        else
                            plotline = plot(t, T, t, V, t, H);
                        end
                        fig.Name = 'energy';
                        Mmin = min([min(V), min(T), min(H)]);
                        Mmax = max([max(V), max(T), max(H)]);
                        ylim([Mmin - 0.1 * abs(Mmax-Mmin), Mmax + 0.1 * abs(Mmax-Mmin)]);
                        title(strcat(integrator_string, ': Energy'));
                        legend('$T$', '$V$', '$H$', 'interpreter', 'latex')
                        xlabel('$t$', 'interpreter', 'latex')
                        ylabel('$\mathrm{[J]}$', 'interpreter', 'latex')

                    case 'angular_momentum'

                        %plots the ang. Mom. about dim-axis over time
                        plotline = plot(t, L(:, :));
                        Lmin = min(L(:));
                        Lmax = max(L(:));
                        ylim([Lmin - 0.1 * abs(Lmax-Lmin), Lmax + 0.1 * abs(Lmax-Lmin)]);
                        fig.Name = 'ang_mom';
                        title(strcat(integrator_string, ': Angular Momentum'));
                        legend('$L_1$', '$L_2$', '$L_3$', 'interpreter', 'latex');
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$L_i(t) \ \mathrm{[J\,s]}$', 'interpreter', 'latex');

                    case 'energy_difference'

                        %plots the increments of Hamiltonian over time
                        plotline = plot(t(1:end-1), diffH, 'color', self.color_scheme{3});
                        max_diff = max(diffH);
                        max_rounded = 10^(real(floor(log10(max_diff))) + 1);
                        min_diff = min(diffH);
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);
                        fig.Name = 'H_diff';
                        title(strcat(integrator_string, ': Energy difference'));
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$H^\mathrm{n+1}-H^\mathrm{n} \ \mathrm{[J]}$', 'interpreter', 'latex');

                    case 'angular_momentum_difference'

                        %plots the increments of ang. Mom. over time about desired axis (dim)
                        plotline = plot(t(1:end-1), diffL(:, :));
                        max_diff = max(max(diffL));
                        max_rounded = 10^(floor(log10(max_diff)) + 1);
                        min_diff = min(min(diffL));
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);
                        fig.Name = 'J_diff';
                        title(strcat(integrator_string, ': Ang. mom. - difference'));
                        xlabel('$t$', 'interpreter', 'latex');
                        legend('$L_1$', '$L_2$', '$L_3$', 'interpreter', 'latex');
                        ylabel('$L_i^\mathrm{n+1}-L_i^\mathrm{n} \ \mathrm{[J\,s]}$', 'interpreter', 'latex');

                    case 'constraint_position'

                        %plots the position constraint and their violations by integration
                        plotline = plot(t, g_pos);
                        max_diff = max(max(g_pos));
                        max_rounded = 10^(floor(log10(max_diff)) + 1);
                        min_diff = min(min(g_pos));
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);

                        fig.Name = 'g_pos';
                        title(strcat(integrator_string, ': Constraint on position level'));
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$g^\mathbf{q}_k(t)$', 'interpreter', 'latex');

                    case 'constraint_velocity'

                        %plots the velocity constraint and their violations by integration
                        plotline = plot(t, g_vel);
                        max_diff = max(max(g_vel));
                        max_rounded = 10^(real(floor(log10(max_diff))+1));
                        min_diff = min(min(g_vel));
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);

                        fig.Name = 'g_vel';
                        title(strcat(integrator_string, ': Constraint on velocity level'));
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$g^\mathbf{v}_k(t)$', 'interpreter', 'latex');

                    otherwise

                        error(strcat('No plotting routine for ', quantity, '  this quantity defined'));

                        end

                        % set linewidth and colorscheme
                        set(plotline, 'linewidth', 1.5);
                        colororder(self.color_scheme);

                end

            end

            %% Function: calculate error for error analysis
                function error = calculate_errors(~, quantity, quantity_ref, num_A, num_B)
                    error = zeros(num_A, num_B);
                    for i = 1:num_A
                        for j = 1:num_B
                            error(i, j) = norm(quantity(:, i, j)-quantity_ref) / norm(quantity_ref);
                        end
                    end
            end

                %% Function: plots convergence-diagram
                    function convergence_plot(self, h_values, y_values, num_A, legend_entries)
                        % plots convergence-diagram y(h) except for last pair of values
                        % since last pair of values might be the reference value
                        % (smallest h)

                        figure()
                        for j = 1:num_A
                            loglog(h_values(1:end-1), y_values(1:end-1, j), '-o', 'Linewidth', 1.2)
                            hold on
                        end
                        xlim([h_values(end-1) * 0.8, h_values(1) * 1.2]);
                        ylim([min(min(y_values(1:end-1, :))) * 0.8, max(max(y_values(1:end-1, :))) * 1.2]);

                        colororder(self.color_scheme);
                        legend(legend_entries)
                        legend('Location', 'southeast')
                        title('relative error')
                        grid on
                        xlabel('h');
                        ylabel('relative error');

                end

                end

            end
