%% Class: Postprocessing
% Distributes animation, computation of postprocessing quantities and plots
classdef Postprocess

    properties

        % Qualitative color palette for all color-vision deficiencies (Okabe and Ito 2008)
        % https://clauswilke.com/dataviz/color-pitfalls.html#not-designing-for-color-vision-deficiency
        color_scheme = {'#E69F00', '#56B4E9', '#CC79A7', '#009E73', '#F0E442', '#0072B2', '#D55E00'};

    end


    methods

        %% Function: initialise the postprocessing
        function self = Postprocess()
            % set default position of figure windows
            set(0, 'defaultfigureposition', [0, 0, 1512, 982])

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
        function this_simulation = compute(~, this_system, this_simulation, this_integrator)

            nDOF = this_system.nDOF;
            m = this_system.mCONSTRAINTS;
            DIM = this_system.DIM;
            d = nDOF / DIM;
            
            NT = int32(this_simulation.T_END/this_simulation.DT) + 1;
            %external_torque = this_system.MASS_MAT;
            %IM = external_torque \ eye(size(external_torque));

            q = this_simulation.z(:, 1:nDOF);
            p = this_simulation.z(:, nDOF+1:2*nDOF);
            v = zeros(NT, nDOF);

            if ismethod(this_system,'get_cartesian_coordinates_center_of_mass')
                x = zeros(NT,this_simulation.DIM);
            end
            
            if this_system.mMixedQuantities > 0
                alpha = this_simulation.z(:,3*nDOF+1:3*nDOF+this_system.mMixedQuantities);
            end

            if m > 0
                lambda = this_simulation.z(:,2*nDOF+1:2*nDOF+m);
            end

            % Check if integration scheme had independent velocity quantities
            if this_simulation.INDI_VELO == true

                % set velocity vector directly
                if strcmp(this_integrator.NAME,'Lagrange_top_ODE')
                    v = this_simulation.z(:, nDOF+1:2*nDOF);
                else
                    v = this_simulation.z(:, 2*nDOF+1:3*nDOF);
                end

                if m > 0 && ~strcmp(this_integrator.NAME,'EML-null')
                    lambda = this_simulation.z(:,3*nDOF+1:3*nDOF+m);
                else 
                    lambda = NaN(NT,m);
                end

                if this_integrator.has_enhanced_constraint_force

                    gamma = this_simulation.z(:,3*nDOF+m+1:3*nDOF+2*m);
                end

            else

                % else: compute by means of momenta
                for j = 1:NT
                    M = this_system.get_mass_matrix(q(j,:));
                    IM =  M \ eye(size(M));
                    v(j, :) = (IM * p(j, :)')';
                end
                

                if this_integrator.has_enhanced_constraint_force
                     gamma = this_simulation.z(:,2*nDOF+m+1:2*nDOF+2*m);
                end

               

            end

            if this_system.mMixedQuantities > 0 && this_simulation.INDI_VELO
                mixed_quantity_difference = zeros(NT,this_system.mMixedQuantities);
                for j = 1:NT
                    mixed_quantity = this_simulation.z(j,3*nDOF+1:3*nDOF+this_system.mMixedQuantities);
                    mixed_quantity_from_position = this_system.mixed_quantity(q(j,:)');
                    mixed_quantity_difference(j) = mixed_quantity - mixed_quantity_from_position;
                end
            end

            % Allocate space for postprocessing quantities
            T = zeros(NT, 1);
            V = zeros(NT, 1);
            H = zeros(NT, 1);
            E = zeros(NT, 1);
            E_kin = zeros(NT, 1);
            D = zeros(NT, 1);
            L = zeros(NT, 3);
            diffH = zeros(NT-1, 1);
            diffE = zeros(NT-1, 1);
            diffL = zeros(NT-1, 3);
            diss_work = zeros(NT-1, 1);
            if DIM == 2
                L = zeros(NT,1);
                diffL = zeros(NT-1,1);
            end
            constraint_position = zeros(NT, m);
            constraint_velocity = zeros(NT, m);
            constraint_forces = zeros(NT-1,nDOF);

                
            external_torque = zeros(NT,3);
            
            alpha_from_q = zeros(NT, 1);

            % Compute quantities for every point in time
            for j = 1:NT
                
                if ismethod(this_system,'get_cartesian_coordinates_center_of_mass')
                    x(j,:) = this_system.get_cartesian_coordinates_center_of_mass(q(j,:));
                else
                    x = [];
                end

                M = this_system.get_mass_matrix(q(j,:));
                
                % Kinetic and potential energy, Hamiltonian
                T(j) = 1 / 2 * v(j, :) * M * v(j, :)'; %in rare cases,
                %compute T with the velocity quantities
                % T(j) = 1/2*p(j,:)*IM*p(j,:)';
                V(j) = this_system.internal_potential(q(j, :)') + this_system.external_potential(q(j, :)');
                
                if this_integrator.compute_potential_from_mixed_quantity
                    V(j) = this_system.internal_potential_from_mixed_quantity(alpha(j, :)') + this_system.external_potential(q(j, :)');
                    alpha_from_q(j) = (q(j, :)*q(j, :)'-1)/2;
                end

                H(j) = T(j) + V(j);
                E(j) = p(j,:)*v(j,:)' - 1 / 2 * v(j, :) * M * v(j, :)' + V(j);
                E_kin(j) = E(j) - V(j);
                D(j) = v(j, :) * this_system.DISS_MAT * v(j, :)';
                F_ext = -this_system.external_potential_gradient(q(j, :)');

                % Constraints on position and velocity level
                if m > 0
                    constraint_position(j, :) = this_system.constraint(q(j, :)')';
                    constraint_velocity(j, :) = (this_system.constraint_gradient(q(j, :)') * v(j, :)')';
                end

                if strcmp(this_simulation.INTEGRATOR, 'GGL_VI_mod')
                    % this specific integration scheme relies upon
                    % secondary constraints with momentum formulation and
                    % has non-aligned velocites
                    IM = M \ eye(size(M));
                    constraint_velocity(j, :) = (this_system.constraint_gradient(q(j, :)') * IM * p(j, :)')';
                end

                % Compute angular momentum
                if DIM == 3
                    if ismethod(this_system,'get_cartesian_angular_momentum_from_momentum')
                       
                        L(j,:) = this_system.get_cartesian_angular_momentum_from_momentum(q(j,:), p(j,:));

                    else
                        r = q(j,:);
                        r_dot_m = p(j,:);

                            for k=1:d
                                L(j,:) = L(j,:)+ cross(r((k - 1)*DIM+1:k*DIM), r_dot_m((k - 1)*DIM+1:k*DIM));
                                %external_torque(j, :) = external_torque(j, :) + cross(q(j, (k - 1)*DIM+1:k*DIM), F_ext((k - 1)*DIM+1:k*DIM));
                            end

                    end
                
                elseif DIM == 2

                    for k =1:d
                        L(j) = L(j) + q(j,(k - 1)*DIM+1)*p(j,(k - 1)*DIM+2) - q(j,(k - 1)*DIM+2)*p(j,(k - 1)*DIM+1);
                    end

                end

            end

            % Compute time-increments in Hamiltonian and angular momentum
            for j = 1:(NT - 1)

                diffH(j) = H(j+1) - H(j);
                diffE(j) = E(j+1) - E(j);
                diffL(j, :) = L(j+1, :) - L(j, :);                
                diss_work(j+1) = diss_work(j)+D(j+1);

                if m > 0
                    if strcmp(this_simulation.INTEGRATOR, 'EMS_std') || strcmp(this_simulation.INTEGRATOR, 'GGL_std') || strcmp(this_simulation.INTEGRATOR, 'MP_ggl') || strcmp(this_simulation.INTEGRATOR, 'MP_std') || strcmp(this_simulation.INTEGRATOR, 'CSE_B') || strcmp(this_simulation.INTEGRATOR, 'EML') || strcmp(this_simulation.INTEGRATOR, 'MP_Livens') || strcmp(this_simulation.INTEGRATOR, 'EML_reduced') || strcmp(this_simulation.INTEGRATOR, 'EML_null')
        
                        constraint_forces(j,:) = (this_system.constraint_gradient(q(j, :)')' * lambda(j+1, :)')';
                        
                    else 
                        t_bar = zeros(nDOF);
                        for l = 1:m
                            t_bar = t_bar + this_system.constraint_hessian(q(j,:), l) * gamma(j,l);
                        end
                        IM = M \ eye(size(M));
                        constraint_forces(j,:) = (this_system.constraint_gradient(q(j, :)')' * lambda(j+1, :)')'  + (t_bar * IM * p(j,:)')' ;
                    end
                end
                
            end

            % Save computed quantities to simulation-object
            this_simulation.x = x;
            this_simulation.T = T;
            this_simulation.V = V;
            this_simulation.H = H;
            this_simulation.E = E;
            this_simulation.E_kin = E_kin;
            this_simulation.D = D;
            this_simulation.diss_work = diss_work;
            this_simulation.L = L;
            this_simulation.Hdiff = diffH;
            this_simulation.Ediff = diffE;
            this_simulation.Ldiff = diffL;
            this_simulation.constraint_position = constraint_position;
            this_simulation.constraint_velocity = constraint_velocity;
            this_simulation.constraint_forces = constraint_forces;
            this_simulation.external_torque = external_torque;
            this_simulation.alpha_from_q = alpha_from_q;
            if this_system.mMixedQuantities > 0 && this_simulation.INDI_VELO
                this_simulation.mixed_quantity_difference = mixed_quantity_difference;
            end
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

                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Exporting results to .mat.                 ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf('     Exporting results to .mat .                 \n');
                fprintf('  \n');

                % Export as .mat file
                save([export_folder, 'results'], 'export_simulation');
                
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Exporting succesful.                 ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
                fprintf('     Exporting succesful.                 \n');
                fprintf('  \n');
                fprintf('**************************************************** \n');

            end

            % Check if plots should be exported as well
            if export_simulation.should_export_figures
                
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Exporting figures.                 ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf('  \n');
                fprintf('     Exporting figures.                 \n');
                fprintf('  \n');

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
                    print(figHandles(i), [export_folder, export_name], '-dpng')

                    % Export to .tikz
                    warning('off')
                    try
                        % add chosen matlab2tikz directory
                        addpath([export_simulation.matlab2tikz_directory,'/src']);
                        % tries to export via matlab2tikz if available
                        
                        matlab2tikz('figurehandle', figHandles(i), 'height', '\figH', 'width', '\figW', 'filename', [export_folder, export_name, '.tikz'], 'showInfo', false, 'floatformat', '%.6g');
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
                        matlab2tikz('figurehandle', figHandles(i), 'height', '\figH', 'width', '\figW', 'filename', [export_folder, export_name, '.tikz'], 'showInfo', false, 'floatformat', '%.6g');
                    end
                    warning('on')

                    %Reset Title and enlarge linewidth
                    set(0, 'currentfigure', figHandles(i));
                    title(titlestring);
                    set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.5);

                end
                
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'     Exporting succesful.                 ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'  ');
                fprintf(export_simulation.log_file_ID, '%s: %s\n', datestr(now, 0),'**************************************************** ');
                fprintf('     Exporting succesful.                 \n');
                fprintf('  \n');
                fprintf('**************************************************** \n');

            end

            % Close log-file
            fclose(export_simulation.log_file_ID);
            % Copy log-file to export folder
            copyfile('metis.log',export_folder)

        end

        %% Function: plot specific quantities
        function plot(self, this_simulation)
            % Function for plotting postprocessing results
            
            x = this_simulation.x;
            H = this_simulation.H;
            E = this_simulation.E;
            E_kin = this_simulation.E_kin;
            V = this_simulation.V;
            T = this_simulation.T;
            t = this_simulation.t;
            L = this_simulation.L;
            diffH = this_simulation.Hdiff;
            diffE = this_simulation.Ediff;
            diffL = this_simulation.Ldiff;
            g_pos = this_simulation.constraint_position;
            g_vel = this_simulation.constraint_velocity;
            C_diff = this_simulation.mixed_quantity_difference;

            % go through desired output quantities
            for i = 1:length(this_simulation.plot_quantities)

                % get current quantity
                fig = figure('Visible', 'off');
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

                    case 'general_energy_function'

                        %plots Energy quantities over time
                        if isnan(E(end))
                            plotline = plot(t(1:end-1), E_kin(1:end-1), t(1:end-1), V(1:end-1), t(1:end-1),  E(1:end-1));
                        else
                            plotline = plot(t, E_kin, t, V, t, E);
                        end
                        fig.Name = 'energy_functional';
                        Mmin = min([min(V), min(E_kin), min(E)]);
                        Mmax = max([max(V), max(E_kin), max(E)]);
                        ylim([Mmin - 0.1 * abs(Mmax-Mmin), Mmax + 0.1 * abs(Mmax-Mmin)]);
                        title(strcat(integrator_string, ': Energy'));
                        legend('$E_{\mathrm{kin}}$', '$V$', '$E$', 'interpreter', 'latex')
                        xlabel('$t$', 'interpreter', 'latex')
                        ylabel('energy', 'interpreter', 'latex')

                    case 'angular_momentum'

                        %plots the ang. Mom. about dim-axis over time
                        plotline = plot(t, L(:, :));
                        Lmin = min(L(:));
                        Lmax = max(L(:));
                        ylim([Lmin - 0.1 * abs(Lmax-Lmin), Lmax + 0.1 * abs(Lmax-Lmin)]);
                        fig.Name = 'ang_mom';
                        title(strcat(integrator_string, ': Angular Momentum'));
                        if this_simulation.DIM == 3
                            legend('$L_1$', '$L_2$', '$L_3$', 'interpreter', 'latex');
                        else
                            legend('$L$', 'interpreter', 'latex');
                        end
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

                    case 'energy_difference_abs'

                        %plots the increments of Hamiltonian over time
                        plotline = plot(t(1:end-1), abs(diffE), 'color', self.color_scheme{3});
                        max_diff = max(diffH);
                        max_rounded = 10^(real(floor(log10(max_diff))) + 1);
                        min_diff = min(diffH);
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);
                        fig.Name = 'E_diff_abs';
                        title(strcat(integrator_string, ': Energy difference'));
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$|E^\mathrm{n+1}-E^\mathrm{n}| \ \mathrm{[J]}$', 'interpreter', 'latex');

                    case 'energy_function_difference'

                        %plots the increments of Hamiltonian over time
                        plotline = plot(t(1:end-1), diffE, 'color', self.color_scheme{3});
                        max_diff = max(diffE);
                        max_rounded = 10^(real(floor(log10(max_diff))) + 1);
                        min_diff = min(diffE);
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);
                        fig.Name = 'E_diff';
                        title(strcat(integrator_string, ': Energy difference'));
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$E^\mathrm{n+1}-E^\mathrm{n}$', 'interpreter', 'latex');

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
                        if this_simulation.DIM == 3
                            legend('$L_1$', '$L_2$', '$L_3$', 'interpreter', 'latex');
                        else
                            legend('$L$', 'interpreter', 'latex');
                        end
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


                    case 'mixed_quantity_difference'

                        %plots the velocity constraint and their violations by integration
                        plotline = plot(t, C_diff);
                        max_diff = max(max(C_diff));
                        max_rounded = 10^(real(floor(log10(max_diff))+1));
                        min_diff = min(min(C_diff));
                        min_rounded = -10^(real(floor(log10(min_diff))) + 1);
                        hold on
                        plot([t(1), t(end)], [max_rounded, max_rounded], 'k--', [t(1), t(end)], [min_rounded, min_rounded], 'k--');
                        ylim([2 * min_rounded, 2 * max_rounded]);

                        fig.Name = 'C_diff';
                        title(strcat(integrator_string, ': mixed quantity difference'));
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$C - C (q) $', 'interpreter', 'latex');


                    case 'cartesian_coordinates_center_of_mass'

                        %plots the ang. Mom. about dim-axis over time
                        plotline = plot(t, x(:, :));
                        Lmin = min(x(:));
                        Lmax = max(x(:));
                        ylim([Lmin - 0.1 * abs(Lmax-Lmin), Lmax + 0.1 * abs(Lmax-Lmin)]);
                        fig.Name = 'cart_coords';
                        title(strcat(integrator_string, ': Cartesian coordinates'));
                        if this_simulation.DIM == 3
                            legend('$x_1$', '$x_2$', '$x_3$', 'interpreter', 'latex');
                        else
                            legend('$x$', 'interpreter', 'latex');
                        end
                        xlabel('$t$', 'interpreter', 'latex');
                        ylabel('$x_i(t) \ \mathrm{[m]}$', 'interpreter', 'latex');


                    otherwise

                        error(strcat('No plotting routine for ', quantity, '  this quantity defined'));

                end

                % set linewidth and colorscheme
                set(plotline, 'linewidth', 1.5);
                colororder(self.color_scheme);

            end

        end

            %% Function: calculate error for error analysis
            function error = calculate_errors(~, this_simulation, quantity, quantity_ref, num_A, num_B)
                    error = zeros(num_A, num_B);

                    for i = 1:num_A
                        for j = 1:num_B
                            if norm(quantity_ref) ~= 0 && ~this_simulation.matrix_error_analysis
                                error(i, j) = norm(quantity(:, i, j)-quantity_ref) / norm(quantity_ref);
                            elseif norm(quantity_ref) == 0
                                error(i, j) = norm(quantity(:, i, j)-quantity_ref);
                            elseif this_simulation.matrix_error_analysis
                                quantity_matrix = reshape(quantity(:,i,j),[3,3]);
                                ref_matrix = reshape(quantity_ref,[3,3]);
                                error(i,j) = norm(ref_matrix*quantity_matrix' - eye(3),'fro');

                            else
                                error('error.')
                            end
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
