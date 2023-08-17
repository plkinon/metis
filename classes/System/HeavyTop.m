%% Class: Gyroscopic Top with steady precession ('heavy top')
%
% A rigid gyroscopic top. Makes use of director formulation, e.g. described
% in [1,2]. Internal constraints plus a constraint which fixes the top to
% the floor. Subject to initial velocities and external acceleration.
%
% References:
% [1]: Betsch, P. and Steinmann, P. Constrained integration of rigid body dynamics.
%      In: Computer Methods in Applied Mechanics and Engineering, 191(3-5): 467–488,
%      2001. doi: 10.1016/S0045-7825(01)00283-3.
%
% [2]: Krenk, S. and Nielsen, M. B. Conservative rigid body dynamics by convected
%      base vectors with implicit constraints. In: Computer Methods in Applied Mechanics
%      and Engineering, 269: 437–453, 2014. doi: 10.1016/j.cma.2013.10.028.

classdef HeavyTop < System
% A rigid gyroscopic top. Makes use of director formulation, e.g. described
% in [1,2]. Internal constraints plus a constraint which fixes the top to
% the floor. Subject to initial velocities and external acceleration.

    %%
    methods

        function self = HeavyTop(CONFIG)

            self.nBODIES = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            % 3 coordinates of center of mass + 3*3 director coordinates
            self.nDOF = 12;
            % ext. acceleration only acts on center of mass
            self.EXT_ACC = [CONFIG.EXT_ACC; zeros(9, 1)];

            % Geometric parameters (correspond to symmetric cone)
            a = 0.1; % total height of the cone
            r = a / 2; % top radius of the cone
            l = 3 * a / 4; % location of center of mass along symmetry axis

            % Principle moments of inertia
            J1 = 3 / 80 * self.MASS * (4 * r^2 + a^2);
            J2 = J1;
            J3 = 3 / 10 * self.MASS * r^2;

            % Entries of convected Euler tensor
            m01 = 1 / 2 * (J2 + J3 - J1);
            m02 = 1 / 2 * (J1 + J3 - J2);
            m03 = 1 / 2 * (J1 + J2 - J3);
            mVec = [self.MASS, self.MASS, self.MASS, m01, m01, m01, m02, m02, m02, m03, m03, m03];

            self.MASS_MAT = diag(mVec);
            self.GEOM = [J1, J2, J3, l];
            % 3 constraints of orthogonality of directors, 3 constraints
            % that directors are normalized to length 1, 3 constraints
            % which fix the cones peak to the origin = 9 in total
            self.mCONSTRAINTS = 9;
            % no internal potential
            self.nPotentialInvariants = 0;
            self.nConstraintInvariants = 9;
            self.nVconstraintInvariants = 9;
            
            self.DISS_MAT = zeros(self.nDOF,self.nDOF);


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
            % given by external acceleration acting on center of mass
            V_ext = -self.MASS * self.EXT_ACC(1:3)' * q(1:3);

        end

        function DV_ext = external_potential_gradient(self, ~)

            DV_ext = zeros(12, 1);
            DV_ext(1:3, 1) = -self.MASS * self.EXT_ACC(1:3);

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
            phi = q(1:self.DIM);
            L = self.GEOM(4);
            d1 = q(self.DIM+1:2*self.DIM);
            d2 = q(2*self.DIM+1:3*self.DIM);
            d3 = q(3*self.DIM+1:4*self.DIM);
            g1 = 0.5 * (d1' * d1 - 1);
            g2 = 0.5 * (d2' * d2 - 1);
            g3 = 0.5 * (d3' * d3 - 1);
            g4 = 0.5 * (d1' * d2);
            g5 = 0.5 * (d1' * d3);
            g6 = 0.5 * (d2' * d3);
            g_cg = phi / L - d3;
            g = [g1; g2; g3; g4; g5; g6; g_cg];

        end

        function Dg = constraint_gradient(self, q)

            % Gradient of constraint w.r.t q
            d1 = q(self.DIM+1:2*self.DIM);
            d2 = q(2*self.DIM+1:3*self.DIM);
            d3 = q(3*self.DIM+1:4*self.DIM);
            null = zeros(1, 3);
            L = self.GEOM(4);
            Dg = 1 / 2 * [null, 2 * d1', null, null; null, null, 2 * d2', null; null, null, null, 2 * d3'; null, d2', d1', null; null, d3', null, d1'; null, null, d3', d2'; 2 * eye(3) / L, zeros(3, 3), zeros(3, 3), -2 * eye(3)];

        end

        function D2g = constraint_hessian(~, ~, m)

            % Hessian of g_1 w.r.t. q
            D2g = zeros(12, 12);
            if m == 1
                D2g(4, 4) = 1;
                D2g(5, 5) = 1;
                D2g(6, 6) = 1;
            elseif m == 2
                D2g(7, 7) = 1;
                D2g(8, 8) = 1;
                D2g(9, 9) = 1;
            elseif m == 3
                D2g(10, 10) = 1;
                D2g(11, 11) = 1;
                D2g(12, 12) = 1;
            elseif m == 4
                D2g(4, 7) = 1 / 2;
                D2g(5, 8) = 1 / 2;
                D2g(6, 9) = 1 / 2;
                D2g(7, 4) = 1 / 2;
                D2g(8, 5) = 1 / 2;
                D2g(9, 6) = 1 / 2;
            elseif m == 5
                D2g(4, 10) = 1 / 2;
                D2g(5, 11) = 1 / 2;
                D2g(6, 12) = 1 / 2;
                D2g(10, 4) = 1 / 2;
                D2g(11, 5) = 1 / 2;
                D2g(12, 6) = 1 / 2;
            elseif m == 6
                D2g(7, 10) = 1 / 2;
                D2g(8, 11) = 1 / 2;
                D2g(9, 12) = 1 / 2;
                D2g(10, 7) = 1 / 2;
                D2g(11, 8) = 1 / 2;
                D2g(12, 9) = 1 / 2;
            end

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

                            L = self.GEOM(4);
                            v = self.MASS_MAT \ p;
                            d1 = q(self.DIM+1:2*self.DIM);
                            d2 = q(2*self.DIM+1:3*self.DIM);
                            d3 = q(3*self.DIM+1:4*self.DIM);
                            v0 = v(1:self.DIM);
                            v1 = v(self.DIM+1:2*self.DIM);
                            v2 = v(2*self.DIM+1:3*self.DIM);
                            v3 = v(3*self.DIM+1:4*self.DIM);

                            if i == 1
                                pi2 = d1' * v1;
                            elseif i == 2
                                pi2 = d2' * v2;
                            elseif i == 3
                                pi2 = d3' * v3;
                            elseif i == 4
                                pi2 = d1' * v2 + d2' * v1;
                            elseif i == 5
                                pi2 = d1' * v3 + d3' * v1;
                            elseif i == 6
                                pi2 = d2' * v3 + d3' * v2;
                            elseif i == 7
                                pi2 = v0(1) / L - v3(1);
                            elseif i == 8
                                pi2 = v0(2) / L - v3(2);
                            elseif i == 9
                                pi2 = v0(3) / L - v3(3);
                            else
                                error('system has only 6 invariants for the constraint.');
                                end

                            end

                            % gradient of the invariant of the velocity constraint w.r.t. q
                                function Dpi2Dq = vConstraint_invariant_gradient_q(self, ~, p, i)
                                    v = self.MASS_MAT \ p;

                                    v1 = v(self.DIM+1:2*self.DIM);
                                    v2 = v(2*self.DIM+1:3*self.DIM);
                                    v3 = v(3*self.DIM+1:4*self.DIM);

                                    if i == 1
                                        Dpi2Dq = [zeros(3, 1); v1; zeros(3, 1); zeros(3, 1)];
                                    elseif i == 2
                                        Dpi2Dq = [zeros(3, 1); zeros(3, 1); v2; zeros(3, 1)];
                                    elseif i == 3
                                        Dpi2Dq = [zeros(3, 1); zeros(3, 1); zeros(3, 1); v3];
                                    elseif i == 4
                                        Dpi2Dq = [zeros(3, 1); v2; v1; zeros(3, 1)];
                                    elseif i == 5
                                        Dpi2Dq = [zeros(3, 1); v3; zeros(3, 1); v1];
                                    elseif i == 6
                                        Dpi2Dq = [zeros(3, 1); zeros(3, 1); v3; v2];
                                    elseif i == 7 || i == 8 || i == 9
                                        Dpi2Dq = zeros(12, 1);
                                    else
                                        error('system has only 6 invariants for the constraint.');
                                        end

                                    end

                                    % gradient of the invariant of the velocity constraint w.r.t. p
                                        function Dpi2Dp = vConstraint_invariant_gradient_p(self, q, ~, i)
                                            M = self.MASS_MAT;
                                            L = self.GEOM(4);
                                            m = M(1, 1);
                                            m1 = M(4, 4);
                                            m2 = M(7, 7);
                                            m3 = M(10, 10);
                                            d1 = q(self.DIM+1:2*self.DIM);
                                            d2 = q(2*self.DIM+1:3*self.DIM);
                                            d3 = q(3*self.DIM+1:4*self.DIM);

                                            if i == 1
                                                Dpi2Dp = [zeros(3, 1); d1 / m1; zeros(3, 1); zeros(3, 1)];
                                            elseif i == 2
                                                Dpi2Dp = [zeros(3, 1); zeros(3, 1); d2 / m2; zeros(3, 1)];
                                            elseif i == 3
                                                Dpi2Dp = [zeros(3, 1); zeros(3, 1); zeros(3, 1); d3 / m3];
                                            elseif i == 4
                                                Dpi2Dp = [zeros(3, 1); d2 / m2; d1 / m1; zeros(3, 1)];
                                            elseif i == 5
                                                Dpi2Dp = [zeros(3, 1); d3 / m3; zeros(3, 1); d1 / m1];
                                            elseif i == 6
                                                Dpi2Dp = [zeros(3, 1); zeros(3, 1); d3 / m3; d2 / m2];
                                            elseif i == 7
                                                Dpi2Dp = [1 / (L * m); 0; 0; zeros(3, 1); zeros(3, 1); -1 / m3; 0; 0];
                                            elseif i == 8
                                                Dpi2Dp = [0; 1 / (L * m); 0; zeros(3, 1); zeros(3, 1); 0; -1 / m3; 0];
                                            elseif i == 9
                                                Dpi2Dp = [0; 0; 1 / (L * m); zeros(3, 1); zeros(3, 1); 0; 0; -1 / m3];
                                            else
                                                error('system has only 6 invariants for the constraint.');
                                                end

                                            end

                                                function D2piDqDp = vConstraint_invariant_hessian_qp(self, ~, ~, i)
                                                    M = self.MASS_MAT;

                                                    m1 = M(4, 4);
                                                    m2 = M(7, 7);
                                                    m3 = M(10, 10);


                                                    if i == 1
                                                        D2piDqDp = [zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), eye(self.DIM) * 1 / m1, zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM)];
                                                    elseif i == 2
                                                        D2piDqDp = [zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), eye(self.DIM) * 1 / m2, zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM)];
                                                    elseif i == 3
                                                        D2piDqDp = [zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), eye(self.DIM) * 1 / m3];
                                                    elseif i == 4
                                                        D2piDqDp = [zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), eye(self.DIM) * 1 / m2, zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), eye(self.DIM) * 1 / m1, zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM)];
                                                    elseif i == 5
                                                        D2piDqDp = [zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), eye(self.DIM) * 1 / m3, zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), eye(self.DIM) * 1 / m1];
                                                    elseif i == 6
                                                        D2piDqDp = [zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), eye(self.DIM) * 1 / m3, zeros(self.DIM); zeros(self.DIM), zeros(self.DIM), zeros(self.DIM), eye(self.DIM) * 1 / m2];
                                                    elseif i == 7 || i == 8 || i == 9
                                                        D2piDqDp = zeros(12, 12);
                                                    end

                                            end

                                                % velocity constraint computed with its invariant
                                                    function gv = Vconstraint_from_invariant(~, pi2, ~)

                                                        gv = pi2;

                                                end

                                                        function DgvDpi = Vconstraint_gradient_from_invariant(~, ~, ~)

                                                            DgvDpi = 1;

                                                    end

                                                        %%%%%%%%%%%%%%%%%%%%%%%%%

                                                            function zeta = constraint_invariant(self, q, i)
                                                                % Constraint on position level
                                                                phi = q(1:self.DIM);
                                                                d1 = q(self.DIM+1:2*self.DIM);
                                                                d2 = q(2*self.DIM+1:3*self.DIM);
                                                                d3 = q(3*self.DIM+1:4*self.DIM);
                                                                L = self.GEOM(4);
                                                                if i == 1
                                                                    zeta = d1' * d1;
                                                                elseif i == 2
                                                                    zeta = d2' * d2;
                                                                elseif i == 3
                                                                    zeta = d3' * d3;
                                                                elseif i == 4
                                                                    zeta = d1' * d2;
                                                                elseif i == 5
                                                                    zeta = d1' * d3;
                                                                elseif i == 6
                                                                    zeta = d2' * d3;
                                                                elseif i == 7
                                                                    zeta = phi(1) / L - d3(1);
                                                                elseif i == 8
                                                                    zeta = phi(2) / L - d3(2);
                                                                elseif i == 9
                                                                    zeta = phi(3) / L - d3(3);
                                                                else
                                                                    error('system has only 6 invariants for the constraint.');
                                                                    end
                                                                end

                                                                    function DzetaDq = constraint_invariant_gradient(self, q, i)
                                                                        d1 = q(self.DIM+1:2*self.DIM);
                                                                        d2 = q(2*self.DIM+1:3*self.DIM);
                                                                        d3 = q(3*self.DIM+1:4*self.DIM);
                                                                        L = self.GEOM(4);

                                                                        if i == 1
                                                                            DzetaDq = [zeros(3, 1); 2 * d1; zeros(3, 1); zeros(3, 1)];
                                                                        elseif i == 2
                                                                            DzetaDq = [zeros(3, 1); zeros(3, 1); 2 * d2; zeros(3, 1)];
                                                                        elseif i == 3
                                                                            DzetaDq = [zeros(3, 1); zeros(3, 1); zeros(3, 1); 2 * d3];
                                                                        elseif i == 4
                                                                            DzetaDq = [zeros(3, 1); d1; d2; zeros(3, 1)];
                                                                        elseif i == 5
                                                                            DzetaDq = [zeros(3, 1); d1; zeros(3, 1); d3];
                                                                        elseif i == 6
                                                                            DzetaDq = [zeros(3, 1); zeros(3, 1); d2; d3];
                                                                        elseif i == 7
                                                                            DzetaDq = [1 / L; 0; 0; zeros(3, 1); zeros(3, 1); -1; 0; 0];
                                                                        elseif i == 8
                                                                            DzetaDq = [0; 1 / L; 0; zeros(3, 1); zeros(3, 1); 0; -1; 0];
                                                                        elseif i == 9
                                                                            DzetaDq = [0; 0; 1 / L; zeros(3, 1); zeros(3, 1); 0; 0; -1];
                                                                        else
                                                                            error('system has only 6 invariants for the constraint.');
                                                                            end
                                                                        end

                                                                        % gradient of the invariant of the position constraint w.r.t. q
                                                                            function D2zetaDq2 = constraint_invariant_hessian(~, ~, i)

                                                                                D2zetaDq2 = zeros(12, 12);
                                                                                if i == 1
                                                                                    D2zetaDq2(4, 4) = 1;
                                                                                    D2zetaDq2(5, 5) = 1;
                                                                                    D2zetaDq2(6, 6) = 1;
                                                                                elseif i == 2
                                                                                    D2zetaDq2(7, 7) = 1;
                                                                                    D2zetaDq2(8, 8) = 1;
                                                                                    D2zetaDq2(9, 9) = 1;
                                                                                elseif i == 3
                                                                                    D2zetaDq2(10, 10) = 1;
                                                                                    D2zetaDq2(11, 11) = 1;
                                                                                    D2zetaDq2(12, 12) = 1;
                                                                                elseif i == 4
                                                                                    D2zetaDq2(4, 7) = 1;
                                                                                    D2zetaDq2(5, 8) = 1;
                                                                                    D2zetaDq2(6, 9) = 1;
                                                                                    D2zetaDq2(7, 4) = 1;
                                                                                    D2zetaDq2(8, 5) = 1;
                                                                                    D2zetaDq2(9, 6) = 1;
                                                                                elseif i == 5
                                                                                    D2zetaDq2(4, 10) = 1;
                                                                                    D2zetaDq2(5, 11) = 1;
                                                                                    D2zetaDq2(6, 12) = 1;
                                                                                    D2zetaDq2(10, 4) = 1;
                                                                                    D2zetaDq2(11, 5) = 1;
                                                                                    D2zetaDq2(12, 6) = 1;
                                                                                elseif i == 6
                                                                                    D2zetaDq2(7, 10) = 1;
                                                                                    D2zetaDq2(8, 11) = 1;
                                                                                    D2zetaDq2(9, 12) = 1;
                                                                                    D2zetaDq2(10, 7) = 1;
                                                                                    D2zetaDq2(11, 8) = 1;
                                                                                    D2zetaDq2(12, 9) = 1;
                                                                                end
                                                                        end

                                                                                function gs = constraint_from_invariant(~, zeta, i)

                                                                                    if i == 1 || i == 2 || i == 3
                                                                                        gs = 0.5 * (zeta - 1);
                                                                                    elseif i == 4 || i == 5 || i == 6
                                                                                        gs = 0.5 * zeta;
                                                                                    elseif i == 7 || i == 8 || i == 9
                                                                                        gs = zeta;
                                                                                    end

                                                                            end

                                                                                    function gs = constraint_gradient_from_invariant(~, ~, i)

                                                                                        gs = 0.5;
                                                                                        if i == 7 || i == 8 || i == 9
                                                                                            gs = 1;
                                                                                        end

                                                                                    end

                                                                                    function analyzed_quantity = hconvergence_set(~, this_simulation)
                                                                                            if strcmp(this_simulation.CONV_QUANTITY,'q')
                                                                                                analyzed_quantity = this_simulation.z(end, 3); %z-coordinate of center of mass
                                                                                            elseif strcmp(this_simulation.CONV_QUANTITY,'p')
                                                                                                analyzed_quantity = this_simulation.z(end, 13:15); %x-coordinate of momentum of center of mass
                                                                                            elseif strcmp(this_simulation.CONV_QUANTITY,'lambda')
                                                                                                analyzed_quantity = this_simulation.z(end, end-15); %LM for external constraint on pos level
                                                                                            
                                                                                            else
                                                                                                error('quantity not yet implement for convergence analysis.')
                                                                                            end

                                                                                        end


                                                                                        function reference_solution = hconvergence_reference(~, this_simulation, analyzed_quantity)

                                                                                            if strcmp(this_simulation.CONV_QUANTITY,'q')    
                                                                                                reference_solution = this_simulation.Q_0(3); %position
                                                                                                
                                                                                            elseif strcmp(this_simulation.CONV_QUANTITY,'p')
                                                                                                %reference_solution = 0; %velocity
                                                                                                reference_solution = analyzed_quantity(:, end, end);
                                                                                            elseif strcmp(this_simulation.CONV_QUANTITY,'lambda')
                                                                                                reference_solution = analyzed_quantity(:, end, end); %LM for external constraint on pos level
                                                                                            else
                                                                                                error('quantity not yet implemented for convergence analysis.')
                                                                                            end

                                                                                             end

                                                                                                function give_animation(self, fig, this_simulation)

                                                                                                    DIM = self.DIM;
                                                                                                    phi = this_simulation.z(:, 1:DIM);
                                                                                                    d1 = this_simulation.z(:, DIM+1:2*DIM);
                                                                                                    d2 = this_simulation.z(:, 2*DIM+1:3*DIM);
                                                                                                    d3 = this_simulation.z(:, 3*DIM+1:4*DIM);
                                                                                                    NT = size(phi, 1);

                                                                                                    axis equal
                                                                                                    minx = min(phi(:, 1)) - 1;
                                                                                                    maxx = max(phi(:, 1)) + 1;
                                                                                                    miny = min(phi(:, 2)) - 1;
                                                                                                    maxy = max(phi(:, 2)) + 1;
                                                                                                    minz = min(phi(:, 3)) - 1;
                                                                                                    maxz = max(phi(:, 3)) + 1;
                                                                                                    axis([minx, maxx, miny, maxy, minz, maxz]);
                                                                                                    xlabel('x');
                                                                                                    ylabel('y');
                                                                                                    zlabel('z');
                                                                                                    grid on;

                                                                                                    xa = phi(1, 1);
                                                                                                    ya = phi(1, 2);
                                                                                                    if DIM == 3
                                                                                                        za = phi(1, 3);
                                                                                                    else
                                                                                                        za = 0;
                                                                                                    end

                                                                                                    for j = 1:NT

                                                                                                        cla(fig);
                                                                                                        hold on

                                                                                                        %% Current position of center of mass
                                                                                                        x0 = phi(j, 1);
                                                                                                        y0 = phi(j, 2);
                                                                                                        z0 = phi(j, 3);
                                                                                                        x1 = x0 + d1(j, 1);
                                                                                                        y1 = y0 + d1(j, 2);
                                                                                                        z1 = z0 + d1(j, 3);
                                                                                                        x2 = x0 + d2(j, 1);
                                                                                                        y2 = y0 + d2(j, 2);
                                                                                                        z2 = z0 + d2(j, 3);
                                                                                                        x3 = x0 + d3(j, 1);
                                                                                                        y3 = y0 + d3(j, 2);
                                                                                                        z3 = z0 + d3(j, 3);

                                                                                                        %% current position of the mass
                                                                                                        hold on
                                                                                                        plot3(x0, y0, z0, 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                                                                                                        plot3(x1, y1, z1, 'mo', 'MarkerSize', 4, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                                                                                                        plot3(x2, y2, z2, 'mo', 'MarkerSize', 4, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                                                                                                        plot3(x3, y3, z3, 'mo', 'MarkerSize', 4, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);

                                                                                                        grid on

                                                                                                        %% current position of the constraint
                                                                                                        px = [0; x0];
                                                                                                        py = [0; y0];
                                                                                                        pz = [0; z0];
                                                                                                        xx = [x0; x1];
                                                                                                        yy = [y0; y1];
                                                                                                        zz = [z0; z1];
                                                                                                        xxx = [x0; x2];
                                                                                                        yyy = [y0; y2];
                                                                                                        zzz = [z0; z2];
                                                                                                        xxxx = [x0; x3];
                                                                                                        yyyy = [y0; y3];
                                                                                                        zzzz = [z0; z3];
                                                                                                        hold on
                                                                                                        plot3(px, py, pz, 'r');
                                                                                                        plot3(xx, yy, zz, 'k', 'linewidth', 1);
                                                                                                        plot3(xxx, yyy, zzz, 'k', 'linewidth', 1);
                                                                                                        plot3(xxxx, yyyy, zzzz, 'k', 'linewidth', 1);

                                                                                                        view(136, 23)

                                                                                                        drawnow

                                                                                                    end

                                                                                            end

                                                                                        end

                                                                                    end
