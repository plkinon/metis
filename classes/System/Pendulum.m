%% Class: Pendulum
%
% A simple pendulum: point mass constrained by a rigid bar to the origin.
% Length of the rod is determined from initial configuration.

classdef Pendulum < System

    %% Pendulum system in 2 or 3 dimensions

    methods

        function self = Pendulum(CONFIG)
            self.nBODIES = 1;
            self.mCONSTRAINTS = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.MASS_MAT = self.MASS * eye(CONFIG.DIM);
            self.nDOF = self.nBODIES * CONFIG.DIM;
            self.mMixedQuantities = 0;
            self.EXT_ACC = repmat(CONFIG.EXT_ACC, self.nBODIES, 1);
            self.GEOM(1) = norm(CONFIG.Q_0(1:CONFIG.DIM)); %length of pendulum
            self.nPotentialInvariants = 0;
            self.nConstraintInvariants = 1;
            self.nVconstraintInvariants = 1;

            self.DISS_MAT = zeros(self.nDOF, self.nDOF);

        end

        function M = get_mass_matrix(self, ~)
            
            M = self.MASS_MAT;

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(~, q, ~)

            Dq_T = zeros(size(q));

        end

        function Dq_T = kinetic_energy_gradient_from_velocity(~, q, ~)

            Dq_T = zeros(size(q));

        end


        function V_ext = external_potential(self, q)
            V_ext = -(self.MASS_MAT * self.EXT_ACC)' * q;
        end

        function DV_ext = external_potential_gradient(self, ~)
            DV_ext = -self.MASS_MAT * self.EXT_ACC;
        end

        function D2V_ext = external_potential_hessian(~, q)
            D2V_ext = zeros(size(q, 1));
        end

        function V_int = internal_potential(~, ~)
            V_int = 0;
        end

        function DV_int = internal_potential_gradient(~, q)
            DV_int = zeros(size(q));
        end

        function D2V_int = internal_potential_hessian(~, q)
            D2V_int = zeros(size(q, 1));
        end


        function g = constraint(self, q)
            % Constraint on position level
            g = 0.5 * (q' * q / self.GEOM(1)^2 - 1);
        end

        function Dg = constraint_gradient(self, q)
            % Gradient of constraint w.r.t q
            Dg = q'/ self.GEOM(1)^2 ;
        end

        function D2g = constraint_hessian(self, q, ~)
            % Hessian of g_1 w.r.t. q
            D2g = eye(size(q, 1))/ self.GEOM(1)^2;
        end

        function [] = potential_invariant(~, ~, ~)

            error('system has only no invariants for the potential.');

            end

                function [] = potential_invariant_gradient(~, ~, ~)

                    error('system has only no invariants for the potential.');

                    end

                    %% Invariant formulations
                    %  e.g. for EMS

                    % invariant of the velocity invariant
                        function pi2 = vConstraint_invariant(self, q, p, i)

                            m = self.MASS;
                            pi2 = q' * p / m;

                    end

                        % gradient of the invariant of the velocity constraint w.r.t. q
                            function Dpi2Dq = vConstraint_invariant_gradient_q(self, ~, p, i)

                                m = self.MASS;
                                Dpi2Dq = p / m;

                        end

                            % gradient of the invariant of the velocity constraint w.r.t. p
                                function Dpi2Dp = vConstraint_invariant_gradient_p(self, q, ~, i)

                                    m = self.MASS;
                                    Dpi2Dp = q / m;

                            end

                                    function D2piDqDp = vConstraint_invariant_hessian_qp(self, ~, ~, i)

                                        m = self.MASS;
                                        D2piDqDp = 1 / m * eye(self.DIM);

                                end

                                    % velocity constraint computed with its invariant
                                        function gv = Vconstraint_from_invariant(self, pi2, ~)

                                            gv = pi2/ self.GEOM(1)^2;

                                    end

                                            function DgvDpi = Vconstraint_gradient_from_invariant(self, ~, ~)
                                                DgvDpi = 1/ self.GEOM(1)^2;
                                        end

                                            %%%%%%%%%%%%%%%%%%%%%%%%%

                                                function zeta = constraint_invariant(~, q, i)

                                                    if i == 1
                                                        zeta = q' * q;
                                                    else
                                                        error('system has only 1 invariant for the constraint.');
                                                        end
                                                    end

                                                        function DzetaDq = constraint_invariant_gradient(~, q, i)

                                                            if i == 1
                                                                DzetaDq = 2 * q';
                                                            else
                                                                error('system has only 1 invariant for the constraint.');
                                                                end
                                                            end

                                                            % gradient of the invariant of the position constraint w.r.t. q
                                                                function D2zetaDq2 = constraint_invariant_hessian(self, ~, i)

                                                                    if i == 1
                                                                        D2zetaDq2 = 2 * eye(self.DIM);
                                                                    else
                                                                        error('system has only 1 invariant for the constraint.');
                                                                        end
                                                                    end

                                                                        function gs = constraint_from_invariant(self, zeta, i)

                                                                            if i == 1
                                                                                gs = 0.5 * (zeta / self.GEOM(1)^2 - 1);
                                                                            else
                                                                                error('system has only 1 invariant for the constraint.');
                                                                                end
                                                                            end

                                                                                function gs = constraint_gradient_from_invariant(self, ~, i)

                                                                                    if i == 1
                                                                                        gs = 0.5/ self.GEOM(1)^2;
                                                                                    else
                                                                                        error('system has only 1 invariant for the constraint.');
                                                                                        end
                                                                                    end

                                                                                        function analyzed_quantity = hconvergence_set(self, this_simulation)

                                                                                            analyzed_quantity = this_simulation.z(end, 1:self.DIM); %position

                                                                                    end


                                                                                            function reference_solution = hconvergence_reference(~, ~, analyzed_quantity)

                                                                                            reference_solution = analyzed_quantity(:, end, end); %position

                                                                                    end

                                                                                            function give_animation(self, fig, this_simulation)

                                                                                                DIM = self.DIM;
                                                                                                q = this_simulation.z(:, 1:DIM);
                                                                                                l = sum(self.GEOM(:));
                                                                                                NT = size(q, 1);

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

                                                                                                for j = 1:NT

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
                                                                                                    plot3(xa, ya, za, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
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
                                                                                                    plot3(x, y, z, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
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

                                                                                end
