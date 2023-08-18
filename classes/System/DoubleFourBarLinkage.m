%% Class: Double Four Bar Linkage
classdef DoubleFourBarLinkage < System
% A system consisting of 5 bars. Makes use of 2D director formulation, 
% e.g. described in [1]. Internal constraints plus a constraint which 
% fixes the top to the floor. Subject to initial velocities and external 
% acceleration. More details can be found in [2].
% References:
% [1]: Betsch, P. and Uhlar, S.: Energy-momentum conserving integration of 
%      multibody dynamics, In: Multibody Syst Dyn, 17: 243–289,
%      2007. doi: 10.1007/s11044-007-9043-9.
%
% [2]: Ivo Roupa, Sérgio B. Gonçalves, Miguel Tavares da Silva: 
%      Kinematics and dynamics of planar multibody systems with fully 
%      Cartesian coordinates and a generic rigid body, 
%      In: Mechanism and Machine Theory, 180: 105134, 2023. 
%      doi: 10.1016/j.mechmachtheory.2022.105134.

    methods

        function self = DoubleFourBarLinkage(CONFIG)
            % This function contructs the desired system

            self.nBODIES = 5; % Number of bodies
            
            self.DIM = CONFIG.DIM; % Number of dimensions
            
            self.MASS = CONFIG.MASS; 
            % 3 coordinates of center of mass + 3*3 director coordinates
            self.nDOF = self.nBODIES*6;
            % ext. acceleration only acts on center of mass
            self.EXT_ACC = [CONFIG.EXT_ACC; zeros(4, 1); CONFIG.EXT_ACC; zeros(4, 1); CONFIG.EXT_ACC; zeros(4, 1); CONFIG.EXT_ACC; zeros(4, 1); CONFIG.EXT_ACC; zeros(4, 1)];

            
            l = 1; % Geometric parameter

            
            J = 1/12*self.MASS*l^2; % Principle moments of inertia

            % Entries of convected Euler tensor
            E1 = J/2;
            E2 = J/2;
            mVec = [self.MASS(1), self.MASS(1), E1(1), E1(1), E2(1), E2(1), ...
                    self.MASS(2), self.MASS(2), E1(2), E1(2), E2(2), E2(2), ...
                    self.MASS(3), self.MASS(3), E1(3), E1(3), E2(3), E2(3), ...
                    self.MASS(4), self.MASS(4), E1(1), E1(4), E2(4), E2(4), ...
                    self.MASS(5), self.MASS(5), E1(5), E1(5), E2(5), E2(5)];

            self.MASS_MAT = diag(mVec);
            self.GEOM = l; 
            self.mCONSTRAINTS = 29;
            % no internal potential
            self.nPotentialInvariants = 0;
            self.nConstraintInvariants = 29;
            self.nVconstraintInvariants = 29;
            
            % dissipation due to damping pots
            eta1 = 0.5;
            eta2 = 0.5;
            Ci = [eye(2) l/2*eye(2) zeros(2,2);
                  zeros(2,2) zeros(2,2) zeros(2,2);
                  l/2*eye(2) l^2/4*eye(2) zeros(2,2)];
            self.DISS_MAT = blkdiag(zeros(6,6), eta1*Ci, zeros(6,6), eta2*Ci, zeros(6,6));


        end

        function M = get_mass_matrix(self, ~)
            % Collects the mass of the system
            
            M = self.MASS_MAT;

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(~, q, ~)
            % Calculates kinetic energy from momentum

            Dq_T = zeros(size(q));

        end

        function Dq_T = kinetic_energy_gradient_from_velocity(~, q, ~)
            % Calculates kinetic energy from velocity

            Dq_T = zeros(size(q));

        end

        function V_ext = external_potential(self, q)
            % given by external acceleration acting on center of mass
            V_ext = 0;
            g = self.EXT_ACC(1:2);

            for i = 1:self.nBODIES
                V_ext = V_ext - self.MASS(i) * g' * q(6*(i-1)+1:6*(i-1)+2);
            end

        end

        function DV_ext = external_potential_gradient(self, ~)
            % calculates the potential gradient given an external force
            g = self.EXT_ACC(1:2);
            DV_ext = zeros(self.nDOF, 1);
            for i = 1:self.nBODIES
                DV_ext(6*(i-1)+1:6*(i-1)+2, 1) = -self.MASS(i) * g;
            end

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

        function F_visc = viscous_forces(self, v)
            
            F_visc = - self.DISS_MAT*v;

        end

        function C = viscous_forces_gradient(self, ~)

            C = - self.DISS_MAT;
        
        end

        function g = constraint(self, q)

            g=[];
            d1  = NaN(self.nBODIES,2);
            phi = NaN(self.nBODIES,2);
            d2  = NaN(self.nBODIES,2);
            % internal rigidity constraints
            for i = 1:self.nBODIES
                phi(i,:) = q(6*(i-1)+1:6*(i-1)+2);
                d1(i,:) = q(6*(i-1)+3:6*(i-1)+4);
                d2(i,:) = q(6*(i-1)+5:6*(i-1)+6);
                g_int1 = 0.5 * (d1(i,:) * d1(i,:)' - 1);
                g_int2 = 0.5 * (d2(i,:) * d2(i,:)' - 1);
                g_int3 = d1(i,:) * d2(i,:)';
                g = [g; g_int1; g_int2; g_int3];
            end
                       
            % external constraints
            l = self.GEOM;
            g_ext1 = phi(1,:) - l/2*d1(1,:); 
            g_ext2 = phi(3,:) + l/2*d1(3,:) - [l 0]; 
            g_ext3 = phi(5,:) + l/2*d1(5,:) - [2*l 0]; 
            g_ext4 = phi(1,:) + l/2*d1(1,:) - (phi(2,:) -l/2*d1(2,:));
            g_ext5 = phi(2,:) + l/2*d1(2,:) - (phi(3,:) -l/2*d1(3,:));
            g_ext6 = phi(2,:) + l/2*d1(2,:) - (phi(4,:) -l/2*d1(4,:));
            g_ext7 = phi(4,:) + l/2*d1(4,:) - (phi(5,:) -l/2*d1(5,:));

            g = [g; g_ext1'; g_ext2'; g_ext3'; g_ext4'; g_ext5'; g_ext6'; g_ext7'];

        end

        function Dg = constraint_gradient(self, q)

            % Gradient of constraint w.r.t q

            d1  = NaN(self.nBODIES,2);
            phi = NaN(self.nBODIES,2);
            d2  = NaN(self.nBODIES,2);
            DG_int_i = cell(1, self.nBODIES);
            % internal rigidity constraints
            for i = 1:self.nBODIES
                phi(i,:) = q(6*(i-1)+1:6*(i-1)+2);
                d1(i,:) = q(6*(i-1)+3:6*(i-1)+4);
                d2(i,:) = q(6*(i-1)+5:6*(i-1)+6);
                DG_int_i{i} = [zeros(1,2) d1(i,:)    zeros(1,2);
                               zeros(1,2) zeros(1,2) d2(i,:);
                               zeros(1,2) d2(i,:)    d1(i,:)];
            end
            Dg_int = blkdiag(DG_int_i{1}, DG_int_i{2}, DG_int_i{3}, DG_int_i{4}, DG_int_i{5});
            
            l = self.GEOM;
            DG_ext_1 = [eye(2) -1/2*l*eye(2) zeros(2,2)];
            DG_ext_left = [eye(2) 1/2*l*eye(2) zeros(2,2)];
            DG_ext_right = [-eye(2) 1/2*l*eye(2) zeros(2,2)];
            
            Dg_ext = [DG_ext_1  zeros(2,6)   zeros(2,6)   zeros(2,6)   zeros(2,6);
                      zeros(2,6)  zeros(2,6)   DG_ext_left   zeros(2,6)   zeros(2,6);
                      zeros(2,6)  zeros(2,6)   zeros(2,6)   zeros(2,6)   DG_ext_left;
                      DG_ext_left DG_ext_right zeros(2,6)   zeros(2,6)   zeros(2,6);
                      zeros(2,6)  DG_ext_left  DG_ext_right zeros(2,6)   zeros(2,6);
                      zeros(2,6)  DG_ext_left  zeros(2,6)   DG_ext_right zeros(2,6);
                      zeros(2,6)  zeros(2,6)   zeros(2,6)   DG_ext_left DG_ext_right];

            Dg = [Dg_int; Dg_ext];


        end

        function D2g = constraint_hessian(~, ~, m)

            % Hessian of g_k w.r.t. q vanish (only linear constraints)

            D2g = zeros(30, 30);
            if ismember(m, [1,4,7,10,13])
                
                D2gint1 = blkdiag(zeros(2,2),eye(2),zeros(2,2));
                if m ==1
                    D2g(1:6,1:6)=D2gint1;
                elseif m==4
                    D2g(7:12,7:12)=D2gint1;
                elseif m==7
                    D2g(13:18,13:18)=D2gint1;
                elseif m==10
                    D2g(19:24,19:24)=D2gint1;
                elseif m==13
                    D2g(25:30,25:30)=D2gint1;
                end

            elseif ismember(m, [2,5,8,11,14])
               
                D2gint2 = blkdiag(zeros(2,2),zeros(2,2), eye(2));
                if m ==2
                    D2g(1:6,1:6)=D2gint2;
                elseif m==5
                    D2g(7:12,7:12)=D2gint2;
                elseif m==8
                    D2g(13:18,13:18)=D2gint2;
                elseif m==11
                    D2g(19:24,19:24)=D2gint2;
                elseif m==14
                    D2g(25:30,25:30)=D2gint2;
                end

            elseif ismember(m, [3,6,9,12,15])
                
                D2gint3 = [zeros(2,2), eye(2), zeros(2,2);
                       eye(2), zeros(2,2), zeros(2,2);
                       zeros(2,2), zeros(2,2), zeros(2,2)];
                if m ==3
                    D2g(1:6,1:6)=D2gint3;
                elseif m==6
                    D2g(7:12,7:12)=D2gint3;
                elseif m==9
                    D2g(13:18,13:18)=D2gint3;
                elseif m==12
                    D2g(19:24,19:24)=D2gint3;
                elseif m==15
                    D2g(25:30,25:30)=D2gint3;
                end

            elseif m > 29
                error('incorrect index of constraint')
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
