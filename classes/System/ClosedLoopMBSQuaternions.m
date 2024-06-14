%% Class: RigidBodyRotatingQuaternions
%
% The rotation of a rigid body in terms of unit-quaternions.
%

classdef ClosedLoopMBSQuaternions < System

    properties

    end

    methods

        function self = ClosedLoopMBSQuaternions(CONFIG)

            self.mCONSTRAINTS = 4*3 + 4*1;
            self.nBODIES = 4;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.nDOF = 4*7;
            self.MASS_MAT = [];
            self.EXT_ACC = CONFIG.EXT_ACC;

           % parameters
            length_bar = 10;
            J_short_axis = self.MASS/12 * (length_bar^2 + 1^2);
            J_long_axis = self.MASS/12 * (1^2 + 1^2);
            self.GEOM = [J_short_axis,J_long_axis,J_short_axis,J_long_axis, J_short_axis, J_short_axis, J_short_axis,J_long_axis,J_short_axis,J_long_axis, J_short_axis, J_short_axis]; % principal inertia [body1, body2, body3, body4]
            self.GEOM = [self.GEOM, length_bar];
            self.DISS_MAT = zeros(self.nDOF,self.nDOF);
            self.nPotentialInvariants = 0;
            self.nKineticInvariants = 8;
            self.nConstraintInvariants = 0;
            self.mMixedQuantities = 0;
            self.isCyclicCoordinate = [];

            %constraint
            self.hasQuadraticConstraints = true;

        end

        function G_q = get_trafo_matrix(~,q)

            if size(q,1) == 1 && size(q,2) == 4
                q = q';
            end

            %extract vector and scalar part form quaternion
            q_vec = q(2:4);
            q_scalar = q(1);

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
            % transformation matrix
            G_q = [-q_vec, q_scalar*eye(3) - q_hat];

        end
        
        function M = get_mass_matrix(self, q)
            
            q_1 = q(4:7);
            q_2 = q(11:14);
            q_3 = q(18:21);
            q_4 = q(25:28);

            G_q1 = self.get_trafo_matrix(q_1);
            G_q2 = self.get_trafo_matrix(q_2);
            G_q3 = self.get_trafo_matrix(q_3);
            G_q4 = self.get_trafo_matrix(q_4);
            
            % classical inertia tensor
            inertia_tensor_1 = diag(self.GEOM(1:3));
            inertia_tensor_2 = diag(self.GEOM(4:6));
            inertia_tensor_3 = diag(self.GEOM(7:9));
            inertia_tensor_4 = diag(self.GEOM(10:12));
            
            % singular mass matrix
            M1 = 4*G_q1'*inertia_tensor_1*G_q1;
            M2 = 4*G_q2'*inertia_tensor_2*G_q2;
            M3 = 4*G_q3'*inertia_tensor_3*G_q3;
            M4 = 4*G_q4'*inertia_tensor_4*G_q4;

            M = blkdiag(self.MASS*eye(3),M1,self.MASS*eye(3),M2,self.MASS*eye(3), M3,self.MASS*eye(3), M4);

        end

        function Dq_T = kinetic_energy_gradient_from_velocity(self, q, v)
            
            q_1 = q(4:7);
            q_2 = q(11:14);
            q_3 = q(18:21);
            q_4 = q(25:28);

            v_1 = v(4:7);
            v_2 = v(11:14);
            v_3 = v(18:21);
            v_4 = v(25:28);

            G_v1 = self.get_trafo_matrix(v_1);
            G_v2 = self.get_trafo_matrix(v_2);
            G_v3 = self.get_trafo_matrix(v_3);
            G_v4 = self.get_trafo_matrix(v_4);
            
            % classical inertia tensor
            inertia_tensor_1 = diag(self.GEOM(1:3));
            inertia_tensor_2 = diag(self.GEOM(4:6));
            inertia_tensor_3 = diag(self.GEOM(7:9));
            inertia_tensor_4 = diag(self.GEOM(10:12));
            
            % singular inertia matrix in v
            M_4_hat_1 = 4*G_v1'*inertia_tensor_1*G_v1;
            M_4_hat_2 = 4*G_v2'*inertia_tensor_2*G_v2;
            M_4_hat_3 = 4*G_v3'*inertia_tensor_3*G_v3;
            M_4_hat_4 = 4*G_v4'*inertia_tensor_4*G_v4;
            
            % partial derivativa of kinetic energy w.r.t. q(uat)
            Dq1_T = M_4_hat_1 * q_1;
            Dq2_T = M_4_hat_2 * q_2;
            Dq3_T = M_4_hat_3 * q_3;
            Dq4_T = M_4_hat_4 * q_4;

            Dq_T = [zeros(3,1); Dq1_T; zeros(3,1); Dq2_T; zeros(3,1); Dq3_T; zeros(3,1); Dq4_T];

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(~, ~, ~)
            
            Dq_T = NaN(4,1);

        end


        function L_single = get_single_angular_momentum_from_momentum(~, q, p)

            if size(q,1) == 1 && size(q,2) == 4
                q = q';
            end

            if size(p,1) == 1 && size(p,2) == 4
                p = p';
            end

            %extract vector and scalar part form quaternion
            q_vec = q(2:4);
            q_scalar = q(1);

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];

            E_q = [-q_vec, q_scalar*eye(3) + q_hat];
            
            L_single = 1/2 * E_q * p ;

        end



        function L = get_cartesian_angular_momentum_from_momentum(self, q, p)

            q_1 = q(4:7);
            q_2 = q(11:14);
            q_3 = q(18:21);
            q_4 = q(25:28);

            p_1 = p(4:7);
            p_2 = p(11:14);
            p_3 = p(18:21);
            p_4 = p(25:28);

            L_single1 = self.get_single_angular_momentum_from_momentum(q_1,p_1);
            L_single2 = self.get_single_angular_momentum_from_momentum(q_2,p_2);
            L_single3 = self.get_single_angular_momentum_from_momentum(q_3,p_3);
            L_single4 = self.get_single_angular_momentum_from_momentum(q_4,p_4);
            
            L_rot = L_single1 + L_single2 + L_single3 + L_single4;

            phi_1 = q(1:3);
            phi_2 = q(8:10);
            phi_3 = q(15:17);
            phi_4 = q(22:24);

            p_phi_1 = p(1:3);
            p_phi_2 = p(8:10);
            p_phi_3 = p(15:17);
            p_phi_4 = p(22:24);

            L_trans = cross(phi_1, p_phi_1) + cross(phi_2, p_phi_2) + cross(phi_3, p_phi_3) + cross(phi_4, p_phi_4);
                        
            L = L_rot + L_trans';

        end

        %% Potential functions

        function V_ext = external_potential(~, ~)
            % External potential
            V_ext = 0;

        end

        function DV_ext = external_potential_gradient(~, q)

            DV_ext = zeros(size(q));

        end

        function D2V_ext = external_potential_hessian(~, q)

            D2V_ext = zeros(size(q),size(q));

        end

        function V_int = internal_potential(~, ~)

            V_int = 0;

        end

        function DV_int = internal_potential_gradient(~, q)

            DV_int = zeros(size(q));

        end

        function D2V_int = internal_potential_hessian(~, q)

            D2V_int = zeros(size(q),size(q));

        end

        function R = get_rotation_matrix(~,q)

            q_vec = q(2:4);
            q_scalar = q(1);

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
            % transformation matrix
            G_q = [-q_vec, q_scalar*eye(3) - q_hat];
            E_q = [-q_vec, q_scalar*eye(3) + q_hat];

            R = E_q*G_q';

        end


        %% Constraint on position level

        function g = constraint(self, q)

            length_bar = self.GEOM(end);

            q_1 = q(4:7);
            q_2 = q(11:14);
            q_3 = q(18:21);
            q_4 = q(25:28);

            phi_1 = q(1:3);
            phi_2 = q(8:10);
            phi_3 = q(15:17);
            phi_4 = q(22:24);
           
            g_int1=1/2*(q_1'*q_1 -1);
            g_int2=1/2*(q_2'*q_2 -1);
            g_int3=1/2*(q_3'*q_3 -1);
            g_int4=1/2*(q_4'*q_4 -1);

            R1 = self.get_rotation_matrix(q_1);
            R2 = self.get_rotation_matrix(q_2);
            R3 = self.get_rotation_matrix(q_3);
            R4 = self.get_rotation_matrix(q_4);

            X1_2 = [0;length_bar/2;0];
            X1_4 = [0;-length_bar/2;0];
            X2_1 = [length_bar/2;0;0];
            X2_3 = [-length_bar/2;0;0];
            X3_2 = [0;length_bar/2;0];
            X3_4 = [0;-length_bar/2;0];
            X4_1 = [length_bar/2;0;0];
            X4_3 = [-length_bar/2;0;0];

            x1_2 = R1 * X1_2;
            x1_4 = R1 * X1_4;
            x2_1 = R2 * X2_1;
            x2_3 = R2 * X2_3;
            x3_2 = R3 * X3_2;
            x3_4 = R3 * X3_4;
            x4_1 = R4 * X4_1;
            x4_3 = R4 * X4_3;

            g_ext12 = phi_1 - phi_2 + x1_2 - x2_1;
            g_ext23 = phi_2 - phi_3 + x2_3 - x3_2;
            g_ext34 = phi_3 - phi_4 + x3_4 - x4_3;
            g_ext41 = phi_4 - phi_1 + x4_1 - x1_4;

            g = [g_int1; g_int2; g_int3; g_int4; g_ext12; g_ext23; g_ext34; g_ext41];

        end

        function X_bar = get_bar_matrix(~,X)

            X_hat = [0, -X(3), X(2);
                    X(3), 0, -X(1);
                    -X(2), X(1), 0];

            X_bar = [0, -X';
                     X, -X_hat];

        end

        function E_q = get_Ematrix(~,q)

            q_vec = q(2:4);
            q_scalar = q(1);

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
            % transformation matrix
            E_q = [-q_vec, q_scalar*eye(3) + q_hat];

        end

        function Dg = constraint_gradient(self, q)

            length_bar = self.GEOM(end);

            q_1 = q(4:7);
            q_2 = q(11:14);
            q_3 = q(18:21);
            q_4 = q(25:28);
           
            Dg_int1=[zeros(1,3), q_1', zeros(1,21)];
            Dg_int2=[zeros(1,10), q_2', zeros(1,14)];
            Dg_int3=[zeros(1,17), q_3', zeros(1,7)];
            Dg_int4=[zeros(1,24), q_4'];

            X1_2 = [0;length_bar/2;0];
            X1_4 = [0;-length_bar/2;0];
            X2_1 = [length_bar/2;0;0];
            X2_3 = [-length_bar/2;0;0];
            X3_2 = [0;length_bar/2;0];
            X3_4 = [0;-length_bar/2;0];
            X4_1 = [length_bar/2;0;0];
            X4_3 = [-length_bar/2;0;0];

            X_bar_12 = self.get_bar_matrix(X1_2);
            X_bar_14 = self.get_bar_matrix(X1_4);
            X_bar_23 = self.get_bar_matrix(X2_3);
            X_bar_21 = self.get_bar_matrix(X2_1);
            X_bar_34 = self.get_bar_matrix(X3_4);
            X_bar_32 = self.get_bar_matrix(X3_2);
            X_bar_43 = self.get_bar_matrix(X4_3);
            X_bar_41 = self.get_bar_matrix(X4_1);

            E_q1 = self.get_Ematrix(q_1);
            E_q2 = self.get_Ematrix(q_2);
            E_q3 = self.get_Ematrix(q_3);
            E_q4 = self.get_Ematrix(q_4);

            Dg_ext12 = [eye(3), 2*E_q1*X_bar_12, -eye(3), -2*E_q2*X_bar_21, zeros(3), zeros(3,4),zeros(3), zeros(3,4)];
            Dg_ext23 = [zeros(3), zeros(3,4), eye(3), 2*E_q2*X_bar_23, -eye(3), -2*E_q3*X_bar_32, zeros(3), zeros(3,4)];
            Dg_ext34 = [zeros(3), zeros(3,4), zeros(3), zeros(3,4), eye(3), 2*E_q3*X_bar_34, -eye(3), -2*E_q4*X_bar_43];
            Dg_ext41 = [-eye(3), -2*E_q1*X_bar_14, zeros(3), zeros(3,4), zeros(3), zeros(3,4), eye(3), 2*E_q4*X_bar_41,];

            Dg = [Dg_int1; Dg_int2; Dg_int3; Dg_int4; Dg_ext12; Dg_ext23; Dg_ext34; Dg_ext41];

        end

        function D2g = constraint_hessian(~, ~, ~)

             D2g=NaN(28,28);

        end

        %% Invariant formulations
        %  e.g. for EMS

        % Invariant of the internal potential
        function pi = potential_invariant(~, ~, ~)

            error('not available.')

        end

            % gradient of potential invariants w.r.t. q
        function DpiDq = potential_invariant_gradient(~, ~, ~)

            error('not available.')
        end

        % gradient of potential invariants w.r.t. q
        function D2piDq = potential_invariant_hessian(~, ~, ~)

           error('not available.')
            
        end

        % internal potential computed with the invariant
        function Vs = potential_from_invariant(~, ~, ~)
            error('not available.')
        end

        % gradient of internal potential w.r.t. the invariant
        function DVsDpi = potential_gradient_from_invariant(~, ~, ~)
          
            error('not available.')
        end

        % Invariant of the internal potential
        function omega = kinetic_energy_invariant(self, q, v, i)

            if i == 1
                
                q_1 = q(4:7);
                G_q1 = self.get_trafo_matrix(q_1);
                v_1 = v(4:7);
                omega = 2*G_q1*v_1;

            elseif i == 2
                
                q_2 = q(11:14);
                G_q2 = self.get_trafo_matrix(q_2);
                v_2 = v(11:14);
                omega = 2*G_q2*v_2;
            
            elseif i == 3
                
                q_3 = q(18:21);
                G_q3 = self.get_trafo_matrix(q_3);
                v_3 = v(18:21);
                omega = 2*G_q3*v_3;

            elseif i == 4

                q_4 = q(25:28);
                G_q4 = self.get_trafo_matrix(q_4);
                v_4 = v(25:28);
                omega = 2*G_q4*v_4;
            
            elseif i == 5

                v_phi_1 = v(1:3);
                omega = 1/2*v_phi_1'*v_phi_1;

            elseif i == 6

                v_phi_2 = v(8:10);
                omega = 1/2*v_phi_2'*v_phi_2;

            elseif i == 7

                v_phi_3 = v(15:17);
                omega = 1/2*v_phi_3'*v_phi_3;

            elseif i == 8

                v_phi_4 = v(22:24);
                omega = 1/2*v_phi_4'*v_phi_4;

            end

        end

            % gradient of potential invariants w.r.t. q
        function DomegaDq = kinetic_energy_invariant_gradient_q(self, ~, v, i)
            
            if i == 1
                
                v_1 = v(4:7);
                G_v1 = self.get_trafo_matrix(v_1);
                DomegaDq = [zeros(3,3), -2*G_v1, zeros(3,21)];

            elseif i == 2
                
                v_2 = v(11:14);
                G_v2 = self.get_trafo_matrix(v_2);
                DomegaDq = [zeros(3,10), -2*G_v2, zeros(3,14)];
            
            elseif i == 3
                
                v_3 = v(18:21);
                G_v3 = self.get_trafo_matrix(v_3);
                DomegaDq = [zeros(3,17), -2*G_v3, zeros(3,7)];

            elseif i == 4

                v_4 = v(25:28);
                G_v4 = self.get_trafo_matrix(v_4);
                DomegaDq = [zeros(3,24), -2*G_v4];
            
            elseif i == 5

                DomegaDq = zeros(1,28);

            elseif i == 6

                DomegaDq = zeros(1,28);

            elseif i == 7

                DomegaDq = zeros(1,28);

            elseif i == 8

                DomegaDq = zeros(1,28);

            end

        end

         % gradient of kinetic energy invariants w.r.t. v
        function DomegaDv = kinetic_energy_invariant_gradient_v(self, q, v, i)
            
            if i == 1
                
                q_1 = q(4:7);
                G_q1 = self.get_trafo_matrix(q_1);
                DomegaDv = [zeros(3,3), 2*G_q1, zeros(3,21)];

            elseif i == 2
                
                q_2 = q(11:14);
                G_q2 = self.get_trafo_matrix(q_2);
                DomegaDv = [zeros(3,10), 2*G_q2, zeros(3,14)];
            
            elseif i == 3
                
                q_3 = q(18:21);
                G_q3 = self.get_trafo_matrix(q_3);
                DomegaDv = [zeros(3,17), 2*G_q3, zeros(3,7)];

            elseif i == 4

                q_4 = q(25:28);
                G_q4 = self.get_trafo_matrix(q_4);
                DomegaDv = [zeros(3,24), 2*G_q4];
            
            elseif i == 5

                v_phi_1 = v(1:3);
                DomegaDv = [v_phi_1', zeros(1,25)];

            elseif i == 6

                v_phi_2 = v(8:10);
                DomegaDv = [zeros(1,7), v_phi_2', zeros(1,18)];

            elseif i == 7

                v_phi_3 = v(15:17);
                DomegaDv = [zeros(1,14), v_phi_3', zeros(1,11)];

            elseif i == 8

                v_phi_4 = v(22:24);
                DomegaDv = [zeros(1,21), v_phi_4', zeros(1,4)];

            end

        end

        % internal potential computed with the invariant
        function Ts = kinetic_energy_from_invariant(self, omega, i)

            
            if i == 1
                
                inertia_tensor_1 = diag(self.GEOM(1:3));
                Ts = 1/2 *omega'*inertia_tensor_1*omega;

            elseif i == 2
                
                inertia_tensor_2 = diag(self.GEOM(4:6));
                Ts = 1/2 *omega'*inertia_tensor_2*omega;
            
            elseif i == 3
                
                inertia_tensor_3 = diag(self.GEOM(7:9));
                Ts = 1/2 *omega'*inertia_tensor_3*omega;

            elseif i == 4

                inertia_tensor_4 = diag(self.GEOM(10:12));
                Ts = 1/2 *omega'*inertia_tensor_4*omega;
            
            elseif i == 5

                Ts = self.MASS*omega;

            elseif i == 6

                Ts = self.MASS*omega;

            elseif i == 7

                Ts = self.MASS*omega;

            elseif i == 8

                Ts = self.MASS*omega;

            end

        end

        % gradient of internal potential w.r.t. the invariant
        function DTsDpi = kinetic_energy_gradient_from_invariant(self, omega, i)
            
            if i == 1
                
                inertia_tensor_1 = diag(self.GEOM(1:3));
                DTsDpi = inertia_tensor_1*omega;

            elseif i == 2
                
                inertia_tensor_2 = diag(self.GEOM(4:6));
                DTsDpi = inertia_tensor_2*omega;
            
            elseif i == 3
                
                inertia_tensor_3 = diag(self.GEOM(7:9));
                DTsDpi = inertia_tensor_3*omega;

            elseif i == 4

                inertia_tensor_4 = diag(self.GEOM(10:12));
                DTsDpi = inertia_tensor_4*omega;
            
            elseif i == 5

                DTsDpi = self.MASS;

            elseif i == 6

                DTsDpi = self.MASS;

            elseif i == 7

                DTsDpi = self.MASS;

            elseif i == 8

                DTsDpi = self.MASS;

            end

        end

        % invariant of the velocity invariant
        function pi2 = vConstraint_invariant(~, ~, ~, ~)

           error('not available.')
        
        end

        % gradient of the invariant of the velocity constraint w.r.t. q
        function Dpi2Dq = vConstraint_invariant_gradient_q(~, ~, ~, ~)

           error('not available.')

    end

    % gradient of the invariant of the velocity constraint w.r.t. p
    function Dpi2Dp = vConstraint_invariant_gradient_p(~, ~, ~, ~)

         error('not available.')

     end

  function D2piDqDp = vConstraint_invariant_hessian_qp(~, ~, ~, ~)

        error('not available.')
  end

% velocity constraint computed with its invariant
    function gv = Vconstraint_from_invariant(~, ~, ~, ~)

      error('not available.')

    end

    function DgvDpi = Vconstraint_gradient_from_invariant(~, ~, ~, ~)

       error('not available.')
    end

    % invariant of the position constraint
    function zeta = constraint_invariant(~, q, ~)

       zeta = q'*q;
       error('not available.')

    end

    % gradient of the invariant of the position constraint w.r.t. q
    function DzetaDq = constraint_invariant_gradient(~, q, ~)

        DzetaDq = 2*q';
        error('not available.')
    end

    % gradient of the invariant of the position constraint w.r.t. q
    function D2zetaDq2 = constraint_invariant_hessian(~, ~, ~)

      D2zetaDq2 = 2*eye(4);
      error('not available.')

    end

    % position constrained computed with its invariant
      function gs = constraint_from_invariant(~, zeta, ~)

           gs = 1/2*(zeta -1);
           error('not available.')
      end

    % gradient of position constrained w.r.t. its invariant
    function Dgs = constraint_gradient_from_invariant(~, ~, ~)

        Dgs = 1/2;
        error('not available.')

    end

    function timefunction = get_timefunction(~,time)

        fm = 100;

        if (time >= 0) && (time <= 0.5)

            timefunction = 2*fm*time;

        elseif (time > 0.5) && (time <= 1)

            timefunction = 2*fm*(1-time);

        elseif time > 1

            timefunction = 0;

        else
            error ("wrong time.")
        end

    end

    function f_ext = get_external_forces(self,q,~,time)
        
        q_1 = q(4:7);
        E_q1 = self.get_Ematrix(q_1);

        f_bar = 8*self.get_timefunction(time)*[1;0;0];
        m_bar = 6*self.get_timefunction(time)*[1;0;0];

        f_ext = [f_bar; 2*E_q1'*m_bar; zeros(21,1)];
        %f_ext = zeros(28,1);
    end

     function analyzed_quantity = hconvergence_set(~, this_simulation)

        if strcmp(this_simulation.CONV_QUANTITY,'q')
            
            error('not available.')

        elseif strcmp(this_simulation.CONV_QUANTITY,'p')
            error('not available.')
        elseif strcmp(this_simulation.CONV_QUANTITY,'lambda')
            error('not available.')
        else
            error('quantity not yet implemented for convergence analysis.')
        end

    end


    function [reference_solution, this_simulation] = hconvergence_reference(~, this_simulation, analyzed_quantity)

        error('not available.')
      
    end

        %% Animation method
            function give_animation(self, fig, this_simulation)
                
                l = self.GEOM(end);
                factor=0.9;
                b=1;
                h=1;

                q = this_simulation.z(:, 1:28);
                NT = size(q, 1);

                q_1 = q(1,4:7);
                q_2 = q(1,11:14);
                q_3 = q(1,18:21);
                q_4 = q(1,25:28);
    
                phi_1 = q(1,1:3);
                phi_2 = q(1,8:10);
                phi_3 = q(1,15:17);
                phi_4 = q(1,22:24);
                
                R1 = self.get_rotation_matrix(q_1');
                R2 = self.get_rotation_matrix(q_2');
                R3 = self.get_rotation_matrix(q_3');
                R4 = self.get_rotation_matrix(q_4');

                P1 = [b/2,factor*l/2,-h/2]';
                P2 = [-b/2,factor*l/2,-h/2]';
                P3 = [-b/2,factor*-l/2,-h/2]';
                P4 = [b/2,factor*-l/2,-h/2]';
                P5 = [b/2,factor*l/2,h/2]';
                P6 = [-b/2,factor*l/2,h/2]';
                P7 = [-b/2,factor*-l/2,h/2]';
                P8 = [b/2,factor*-l/2,h/2]';

                Q1 = [factor*l/2,b/2,-h/2]';
                Q2 = [factor*-l/2,b/2,-h/2]';
                Q3 = [factor*-l/2,-b/2,-h/2]';
                Q4 = [factor*l/2,-b/2,-h/2]';
                Q5 = [factor*l/2,b/2,h/2]';
                Q6 = [factor*-l/2,b/2,h/2]';
                Q7 = [factor*-l/2,-b/2,h/2]';
                Q8 = [factor*l/2,-b/2,h/2]';

                X12 = [0,l/2,0];
                X14 = [0,-l/2,0];
                X32 = [0,l/2,0];
                X34 = [0,-l/2,0];

                x12 = R1*X12';
                x14 = R1*X14';
                x32 = R3*X32';
                x34 = R3*X34';

                p1a = R1*P1;
                p2a = R1*P2;
                p3a = R1*P3;
                p4a = R1*P4;
                p5a = R1*P5;
                p6a = R1*P6;
                p7a = R1*P7;
                p8a = R1*P8;

                q1a = R2*Q1;
                q2a = R2*Q2;
                q3a = R2*Q3;
                q4a = R2*Q4;
                q5a = R2*Q5;
                q6a = R2*Q6;
                q7a = R2*Q7;
                q8a = R2*Q8;

                r1a = R3*P1;
                r2a = R3*P2;
                r3a = R3*P3;
                r4a = R3*P4;
                r5a = R3*P5;
                r6a = R3*P6;
                r7a = R3*P7;
                r8a = R3*P8;

                s1a = R4*Q1;
                s2a = R4*Q2;
                s3a = R4*Q3;
                s4a = R4*Q4;
                s5a = R4*Q5;
                s6a = R4*Q6;
                s7a = R4*Q7;
                s8a = R4*Q8;

                xa = phi_1(1,1) + [p1a(1);p2a(1);p3a(1);p4a(1);p1a(1);p5a(1);p6a(1);p7a(1);p8a(1);p5a(1)];
                ya = phi_1(1,2) + [p1a(2);p2a(2);p3a(2);p4a(2);p1a(2);p5a(2);p6a(2);p7a(2);p8a(2);p5a(2)];
                za = phi_1(1,3) + [p1a(3);p2a(3);p3a(3);p4a(3);p1a(3);p5a(3);p6a(3);p7a(3);p8a(3);p5a(3)];

                x2a = phi_1(1,1) + [p4a(1);p8a(1);p7a(1);p3a(1);p2a(1);p6a(1)];
                y2a = phi_1(1,2) + [p4a(2);p8a(2);p7a(2);p3a(2);p2a(2);p6a(2)];
                z2a = phi_1(1,3) + [p4a(3);p8a(3);p7a(3);p3a(3);p2a(3);p6a(3)];

                x3a = phi_2(1,1) + [q1a(1);q2a(1);q3a(1);q4a(1);q1a(1);q5a(1);q6a(1);q7a(1);q8a(1);q5a(1)];
                y3a = phi_2(1,2) + [q1a(2);q2a(2);q3a(2);q4a(2);q1a(2);q5a(2);q6a(2);q7a(2);q8a(2);q5a(2)];
                z3a = phi_2(1,3) + [q1a(3);q2a(3);q3a(3);q4a(3);q1a(3);q5a(3);q6a(3);q7a(3);q8a(3);q5a(3)];

                x4a = phi_2(1,1) + [q4a(1);q8a(1);q7a(1);q3a(1);q2a(1);q6a(1)];
                y4a = phi_2(1,2) + [q4a(2);q8a(2);q7a(2);q3a(2);q2a(2);q6a(2)];
                z4a = phi_2(1,3) + [q4a(3);q8a(3);q7a(3);q3a(3);q2a(3);q6a(3)];

                x5a = phi_3(1,1) + [r1a(1);r2a(1);r3a(1);r4a(1);r1a(1);r5a(1);r6a(1);r7a(1);r8a(1);r5a(1)];
                y5a = phi_3(1,2) + [r1a(2);r2a(2);r3a(2);r4a(2);r1a(2);r5a(2);r6a(2);r7a(2);r8a(2);r5a(2)];
                z5a = phi_3(1,3) + [r1a(3);r2a(3);r3a(3);r4a(3);r1a(3);r5a(3);r6a(3);r7a(3);r8a(3);r5a(3)];

                x6a = phi_3(1,1) + [r4a(1);r8a(1);r7a(1);r3a(1);r2a(1);r6a(1)];
                y6a = phi_3(1,2) + [r4a(2);r8a(2);r7a(2);r3a(2);r2a(2);r6a(2)];
                z6a = phi_3(1,3) + [r4a(3);r8a(3);r7a(3);r3a(3);r2a(3);r6a(3)];

                x7a = phi_4(1,1) + [s1a(1);s2a(1);s3a(1);s4a(1);s1a(1);s5a(1);s6a(1);s7a(1);s8a(1);s5a(1)];
                y7a = phi_4(1,2) + [s1a(2);s2a(2);s3a(2);s4a(2);s1a(2);s5a(2);s6a(2);s7a(2);s8a(2);s5a(2)];
                z7a = phi_4(1,3) + [s1a(3);s2a(3);s3a(3);s4a(3);s1a(3);s5a(3);s6a(3);s7a(3);s8a(3);s5a(3)];

                x8a = phi_4(1,1) + [s4a(1);s8a(1);s7a(1);s3a(1);s2a(1);s6a(1)];
                y8a = phi_4(1,2) + [s4a(2);s8a(2);s7a(2);s3a(2);s2a(2);s6a(2)];
                z8a = phi_4(1,3) + [s4a(3);s8a(3);s7a(3);s3a(3);s2a(3);s6a(3)];

                new_x_1 = phi_1(1,1) + [p1a(1);p2a(1);p3a(1);p4a(1)];
                new_y_1 = phi_1(1,2) + [p1a(2);p2a(2);p3a(2);p4a(2)];
                new_z_1 = phi_1(1,3) + [p1a(3);p2a(3);p3a(3);p4a(3)];

                new_x_2 = phi_1(1,1) + [p1a(1);p2a(1);p6a(1);p5a(1)];
                new_y_2 = phi_1(1,2) + [p1a(2);p2a(2);p6a(2);p5a(2)];
                new_z_2 = phi_1(1,3) + [p1a(3);p2a(3);p6a(3);p5a(3)];

                new_x_3 = phi_1(1,1) + [p5a(1);p6a(1);p7a(1);p8a(1)];
                new_y_3 = phi_1(1,2) + [p5a(2);p6a(2);p7a(2);p8a(2)];
                new_z_3 = phi_1(1,3) + [p5a(3);p6a(3);p7a(3);p8a(3)];

                new_x_4 = phi_1(1,1) + [p4a(1);p3a(1);p7a(1);p8a(1)];
                new_y_4 = phi_1(1,2) + [p4a(2);p3a(2);p7a(2);p8a(2)];
                new_z_4 = phi_1(1,3) + [p4a(3);p3a(3);p7a(3);p8a(3)];

                new_x_5 = phi_1(1,1) + [p5a(1);p8a(1);p4a(1);p1a(1)];
                new_y_5 = phi_1(1,2) + [p5a(2);p8a(2);p4a(2);p1a(2)];
                new_z_5 = phi_1(1,3) + [p5a(3);p8a(3);p4a(3);p1a(3)];

                new_x_6 = phi_1(1,1) + [p2a(1);p3a(1);p7a(1);p6a(1)];
                new_y_6 = phi_1(1,2) + [p2a(2);p3a(2);p7a(2);p6a(2)];
                new_z_6 = phi_1(1,3) + [p2a(3);p3a(3);p7a(3);p6a(3)];

                new_x_7 = phi_2(1,1) + [q1a(1);q2a(1);q3a(1);q4a(1)];
                new_y_7 = phi_2(1,2) + [q1a(2);q2a(2);q3a(2);q4a(2)];
                new_z_7 = phi_2(1,3) + [q1a(3);q2a(3);q3a(3);q4a(3)];

                new_x_8 = phi_2(1,1) + [q1a(1);q2a(1);q6a(1);q5a(1)];
                new_y_8 = phi_2(1,2) + [q1a(2);q2a(2);q6a(2);q5a(2)];
                new_z_8 = phi_2(1,3) + [q1a(3);q2a(3);q6a(3);q5a(3)];

                new_x_9 = phi_2(1,1) + [q5a(1);q6a(1);q7a(1);q8a(1)];
                new_y_9 = phi_2(1,2) + [q5a(2);q6a(2);q7a(2);q8a(2)];
                new_z_9 = phi_2(1,3) + [q5a(3);q6a(3);q7a(3);q8a(3)];

                new_x_10 = phi_2(1,1) + [q4a(1);q3a(1);q7a(1);q8a(1)];
                new_y_10 = phi_2(1,2) + [q4a(2);q3a(2);q7a(2);q8a(2)];
                new_z_10 = phi_2(1,3) + [q4a(3);q3a(3);q7a(3);q8a(3)];

                new_x_11 = phi_2(1,1) + [q5a(1);q8a(1);q4a(1);q1a(1)];
                new_y_11 = phi_2(1,2) + [q5a(2);q8a(2);q4a(2);q1a(2)];
                new_z_11 = phi_2(1,3) + [q5a(3);q8a(3);q4a(3);q1a(3)];

                new_x_12 = phi_2(1,1) + [q2a(1);q3a(1);q7a(1);q6a(1)];
                new_y_12 = phi_2(1,2) + [q2a(2);q3a(2);q7a(2);q6a(2)];
                new_z_12 = phi_2(1,3) + [q2a(3);q3a(3);q7a(3);q6a(3)];

                new_x_13 = phi_3(1,1) + [r1a(1);r2a(1);r3a(1);r4a(1)];
                new_y_13 = phi_3(1,2) + [r1a(2);r2a(2);r3a(2);r4a(2)];
                new_z_13 = phi_3(1,3) + [r1a(3);r2a(3);r3a(3);r4a(3)];

                new_x_14 = phi_3(1,1) + [r1a(1);r2a(1);r6a(1);r5a(1)];
                new_y_14 = phi_3(1,2) + [r1a(2);r2a(2);r6a(2);r5a(2)];
                new_z_14 = phi_3(1,3) + [r1a(3);r2a(3);r6a(3);r5a(3)];

                new_x_15 = phi_3(1,1) + [r5a(1);r6a(1);r7a(1);r8a(1)];
                new_y_15 = phi_3(1,2) + [r5a(2);r6a(2);r7a(2);r8a(2)];
                new_z_15 = phi_3(1,3) + [r5a(3);r6a(3);r7a(3);r8a(3)];

                new_x_16 = phi_3(1,1) + [r4a(1);r3a(1);r7a(1);r8a(1)];
                new_y_16 = phi_3(1,2) + [r4a(2);r3a(2);r7a(2);r8a(2)];
                new_z_16 = phi_3(1,3) + [r4a(3);r3a(3);r7a(3);r8a(3)];

                new_x_17 = phi_3(1,1) + [r5a(1);r8a(1);r4a(1);r1a(1)];
                new_y_17 = phi_3(1,2) + [r5a(2);r8a(2);r4a(2);r1a(2)];
                new_z_17 = phi_3(1,3) + [r5a(3);r8a(3);r4a(3);r1a(3)];

                new_x_18 = phi_3(1,1) + [r2a(1);r3a(1);r7a(1);r6a(1)];
                new_y_18 = phi_3(1,2) + [r2a(2);r3a(2);r7a(2);r6a(2)];
                new_z_18 = phi_3(1,3) + [r2a(3);r3a(3);r7a(3);r6a(3)];

                new_x_19 = phi_4(1,1) + [s1a(1);s2a(1);s3a(1);s4a(1)];
                new_y_19 = phi_4(1,2) + [s1a(2);s2a(2);s3a(2);s4a(2)];
                new_z_19 = phi_4(1,3) + [s1a(3);s2a(3);s3a(3);s4a(3)];

                new_x_20 = phi_4(1,1) + [s1a(1);s2a(1);s6a(1);s5a(1)];
                new_y_20 = phi_4(1,2) + [s1a(2);s2a(2);s6a(2);s5a(2)];
                new_z_20 = phi_4(1,3) + [s1a(3);s2a(3);s6a(3);s5a(3)];

                new_x_21 = phi_4(1,1) + [s5a(1);s6a(1);s7a(1);s8a(1)];
                new_y_21 = phi_4(1,2) + [s5a(2);s6a(2);s7a(2);s8a(2)];
                new_z_21 = phi_4(1,3) + [s5a(3);s6a(3);s7a(3);s8a(3)];

                new_x_22 = phi_4(1,1) + [s4a(1);s3a(1);s7a(1);s8a(1)];
                new_y_22 = phi_4(1,2) + [s4a(2);s3a(2);s7a(2);s8a(2)];
                new_z_22 = phi_4(1,3) + [s4a(3);s3a(3);s7a(3);s8a(3)];

                new_x_23 = phi_4(1,1) + [s5a(1);s8a(1);s4a(1);s1a(1)];
                new_y_23 = phi_4(1,2) + [s5a(2);s8a(2);s4a(2);s1a(2)];
                new_z_23 = phi_4(1,3) + [s5a(3);s8a(3);s4a(3);s1a(3)];

                new_x_24 = phi_4(1,1) + [s2a(1);s3a(1);s7a(1);s6a(1)];
                new_y_24 = phi_4(1,2) + [s2a(2);s3a(2);s7a(2);s6a(2)];
                new_z_24 = phi_4(1,3) + [s2a(3);s3a(3);s7a(3);s6a(3)];

                plot3(xa, ya, za, 'k', LineWidth=1.2);
                hold on
                plot3(x2a, y2a, z2a, 'k', LineWidth=1.2);
                plot3(x3a, y3a, z3a, 'k', LineWidth=1.2);
                plot3(x4a, y4a, z4a, 'k', LineWidth=1.2);
                plot3(x5a, y5a, z5a, 'k', LineWidth=1.2);
                plot3(x6a, y6a, z6a, 'k', LineWidth=1.2);
                plot3(x7a, y7a, z7a, 'k', LineWidth=1.2);
                plot3(x8a, y8a, z8a, 'k', LineWidth=1.2);

                grayColor = [.7 .7 .7];

                fill3(new_x_1,new_y_1, new_z_1,"white")
                fill3(new_x_2,new_y_2, new_z_2,"white")
                fill3(new_x_3,new_y_3, new_z_3,"white")
                fill3(new_x_4,new_y_4, new_z_4,"white")
                fill3(new_x_5,new_y_5, new_z_5,"white")
                fill3(new_x_6,new_y_6, new_z_6,"white")
                fill3(new_x_7,new_y_7, new_z_7,"white")
                fill3(new_x_8,new_y_8, new_z_8,"white")
                fill3(new_x_9,new_y_9, new_z_9,"white")
                fill3(new_x_10,new_y_10, new_z_10,"white")
                fill3(new_x_11,new_y_11, new_z_11,"white")
                fill3(new_x_12,new_y_12, new_z_12,"white")
                fill3(new_x_13,new_y_13, new_z_13,"white")
                fill3(new_x_14,new_y_14, new_z_14,"white")
                fill3(new_x_15,new_y_15, new_z_15,"white")
                fill3(new_x_16,new_y_16, new_z_16,"white")
                fill3(new_x_17,new_y_17, new_z_17,"white")
                fill3(new_x_18,new_y_18, new_z_18,"white")
                fill3(new_x_19,new_y_19, new_z_19,"white")
                fill3(new_x_20,new_y_20, new_z_20,"white")
                fill3(new_x_21,new_y_21, new_z_21,"white")
                fill3(new_x_22,new_y_22, new_z_22,"white")
                fill3(new_x_23,new_y_23, new_z_23,"white")
                fill3(new_x_24,new_y_24, new_z_24,"white")
               
                plot3(phi_1(1,1) + x12(1), phi_1(1,2) + x12(2), phi_1(1,3)+ x12(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);
                plot3(phi_1(1,1) + x14(1), phi_1(1,2) + x14(2), phi_1(1,3)+ x14(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);
                plot3(phi_3(1,1) + x32(1), phi_3(1,2) + x32(2), phi_3(1,3)+ x32(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);
                plot3(phi_3(1,1) + x34(1), phi_3(1,2) + x34(2), phi_3(1,3)+ x34(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);

                axis equal
                % xmin = -l - 1;
                % xmax = l + 1;
                % ymin = -l -1;
                % ymax = l + 1;
                % zmin =-l - 1;
                % zmax =l +1;
                % axis([xmin, xmax, ymin, ymax, zmin, zmax]);

                grid on;
                axis off

                for j = 1:(NT-1)/5:NT

                    %cla(fig);
                    hold on

                    %% Current position
                    q_1_current = q(j,4:7);
                    q_2_current = q(j,11:14);
                    q_3_current = q(j,18:21);
                    q_4_current = q(j,25:28);
        
                    phi_1_current = q(j,1:3);
                    phi_2_current = q(j,8:10);
                    phi_3_current = q(j,15:17);
                    phi_4_current = q(j,22:24);
                    
                    R1_current = self.get_rotation_matrix(q_1_current');
                    R2_current = self.get_rotation_matrix(q_2_current');
                    R3_current = self.get_rotation_matrix(q_3_current');
                    R4_current = self.get_rotation_matrix(q_4_current');
                        
                    p1a = R1_current*P1;
                    p2a = R1_current*P2;
                    p3a = R1_current*P3;
                    p4a = R1_current*P4;
                    p5a = R1_current*P5;
                    p6a = R1_current*P6;
                    p7a = R1_current*P7;
                    p8a = R1_current*P8;

                    q1a = R2_current*Q1;
                    q2a = R2_current*Q2;
                    q3a = R2_current*Q3;
                    q4a = R2_current*Q4;
                    q5a = R2_current*Q5;
                    q6a = R2_current*Q6;
                    q7a = R2_current*Q7;
                    q8a = R2_current*Q8;

                    r1a = R3_current*P1;
                    r2a = R3_current*P2;
                    r3a = R3_current*P3;
                    r4a = R3_current*P4;
                    r5a = R3_current*P5;
                    r6a = R3_current*P6;
                    r7a = R3_current*P7;
                    r8a = R3_current*P8;

                    s1a = R4_current*Q1;
                    s2a = R4_current*Q2;
                    s3a = R4_current*Q3;
                    s4a = R4_current*Q4;
                    s5a = R4_current*Q5;
                    s6a = R4_current*Q6;
                    s7a = R4_current*Q7;
                    s8a = R4_current*Q8;

                    x12 = R1_current*X12';
                    x14 = R1_current*X14';
                    x32 = R3_current*X32';
                    x34 = R3_current*X34';
   
                    x = phi_1_current(1,1) + [p1a(1);p2a(1);p3a(1);p4a(1);p1a(1);p5a(1);p6a(1);p7a(1);p8a(1);p5a(1)];
                    y = phi_1_current(1,2) + [p1a(2);p2a(2);p3a(2);p4a(2);p1a(2);p5a(2);p6a(2);p7a(2);p8a(2);p5a(2)];
                    z = phi_1_current(1,3) + [p1a(3);p2a(3);p3a(3);p4a(3);p1a(3);p5a(3);p6a(3);p7a(3);p8a(3);p5a(3)];
                    
                    x2 = phi_1_current(1,1) + [p4a(1);p8a(1);p7a(1);p3a(1);p2a(1);p6a(1)];
                    y2 = phi_1_current(1,2) + [p4a(2);p8a(2);p7a(2);p3a(2);p2a(2);p6a(2)];
                    z2 = phi_1_current(1,3) + [p4a(3);p8a(3);p7a(3);p3a(3);p2a(3);p6a(3)];

                    x3 = phi_2_current(1,1) + [q1a(1);q2a(1);q3a(1);q4a(1);q1a(1);q5a(1);q6a(1);q7a(1);q8a(1);q5a(1)];
                    y3 = phi_2_current(1,2) + [q1a(2);q2a(2);q3a(2);q4a(2);q1a(2);q5a(2);q6a(2);q7a(2);q8a(2);q5a(2)];
                    z3 = phi_2_current(1,3) + [q1a(3);q2a(3);q3a(3);q4a(3);q1a(3);q5a(3);q6a(3);q7a(3);q8a(3);q5a(3)];
    
                    x4 = phi_2_current(1,1) + [q4a(1);q8a(1);q7a(1);q3a(1);q2a(1);q6a(1)];
                    y4 = phi_2_current(1,2) + [q4a(2);q8a(2);q7a(2);q3a(2);q2a(2);q6a(2)];
                    z4 = phi_2_current(1,3) + [q4a(3);q8a(3);q7a(3);q3a(3);q2a(3);q6a(3)];

                    x5 = phi_3_current(1,1) + [r1a(1);r2a(1);r3a(1);r4a(1);r1a(1);r5a(1);r6a(1);r7a(1);r8a(1);r5a(1)];
                    y5 = phi_3_current(1,2) + [r1a(2);r2a(2);r3a(2);r4a(2);r1a(2);r5a(2);r6a(2);r7a(2);r8a(2);r5a(2)];
                    z5 = phi_3_current(1,3) + [r1a(3);r2a(3);r3a(3);r4a(3);r1a(3);r5a(3);r6a(3);r7a(3);r8a(3);r5a(3)];
    
                    x6 = phi_3_current(1,1) + [r4a(1);r8a(1);r7a(1);r3a(1);r2a(1);r6a(1)];
                    y6 = phi_3_current(1,2) + [r4a(2);r8a(2);r7a(2);r3a(2);r2a(2);r6a(2)];
                    z6 = phi_3_current(1,3) + [r4a(3);r8a(3);r7a(3);r3a(3);r2a(3);r6a(3)];
    
                    x7 = phi_4_current(1,1) + [s1a(1);s2a(1);s3a(1);s4a(1);s1a(1);s5a(1);s6a(1);s7a(1);s8a(1);s5a(1)];
                    y7 = phi_4_current(1,2) + [s1a(2);s2a(2);s3a(2);s4a(2);s1a(2);s5a(2);s6a(2);s7a(2);s8a(2);s5a(2)];
                    z7 = phi_4_current(1,3) + [s1a(3);s2a(3);s3a(3);s4a(3);s1a(3);s5a(3);s6a(3);s7a(3);s8a(3);s5a(3)];
    
                    x8 = phi_4_current(1,1) + [s4a(1);s8a(1);s7a(1);s3a(1);s2a(1);s6a(1)];
                    y8 = phi_4_current(1,2) + [s4a(2);s8a(2);s7a(2);s3a(2);s2a(2);s6a(2)];
                    z8 = phi_4_current(1,3) + [s4a(3);s8a(3);s7a(3);s3a(3);s2a(3);s6a(3)];
                    
                    new_x_1 = phi_1_current(1,1) + [p1a(1);p2a(1);p3a(1);p4a(1)];
                    new_y_1 = phi_1_current(1,2) + [p1a(2);p2a(2);p3a(2);p4a(2)];
                    new_z_1 = phi_1_current(1,3) + [p1a(3);p2a(3);p3a(3);p4a(3)];
    
                    new_x_2 = phi_1_current(1,1) + [p1a(1);p2a(1);p6a(1);p5a(1)];
                    new_y_2 = phi_1_current(1,2) + [p1a(2);p2a(2);p6a(2);p5a(2)];
                    new_z_2 = phi_1_current(1,3) + [p1a(3);p2a(3);p6a(3);p5a(3)];
    
                    new_x_3 = phi_1_current(1,1) + [p5a(1);p6a(1);p7a(1);p8a(1)];
                    new_y_3 = phi_1_current(1,2) + [p5a(2);p6a(2);p7a(2);p8a(2)];
                    new_z_3 = phi_1_current(1,3) + [p5a(3);p6a(3);p7a(3);p8a(3)];
    
                    new_x_4 = phi_1_current(1,1) + [p4a(1);p3a(1);p7a(1);p8a(1)];
                    new_y_4 = phi_1_current(1,2) + [p4a(2);p3a(2);p7a(2);p8a(2)];
                    new_z_4 = phi_1_current(1,3) + [p4a(3);p3a(3);p7a(3);p8a(3)];
    
                    new_x_5 = phi_1_current(1,1) + [p5a(1);p8a(1);p4a(1);p1a(1)];
                    new_y_5 = phi_1_current(1,2) + [p5a(2);p8a(2);p4a(2);p1a(2)];
                    new_z_5 = phi_1_current(1,3) + [p5a(3);p8a(3);p4a(3);p1a(3)];
    
                    new_x_6 = phi_1_current(1,1) + [p2a(1);p3a(1);p7a(1);p6a(1)];
                    new_y_6 = phi_1_current(1,2) + [p2a(2);p3a(2);p7a(2);p6a(2)];
                    new_z_6 = phi_1_current(1,3) + [p2a(3);p3a(3);p7a(3);p6a(3)];

                    new_x_7 = phi_2_current(1,1) + [q1a(1);q2a(1);q3a(1);q4a(1)];
                    new_y_7 = phi_2_current(1,2) + [q1a(2);q2a(2);q3a(2);q4a(2)];
                    new_z_7 = phi_2_current(1,3) + [q1a(3);q2a(3);q3a(3);q4a(3)];
    
                    new_x_8 = phi_2_current(1,1) + [q1a(1);q2a(1);q6a(1);q5a(1)];
                    new_y_8 = phi_2_current(1,2) + [q1a(2);q2a(2);q6a(2);q5a(2)];
                    new_z_8 = phi_2_current(1,3) + [q1a(3);q2a(3);q6a(3);q5a(3)];
    
                    new_x_9 = phi_2_current(1,1) + [q5a(1);q6a(1);q7a(1);q8a(1)];
                    new_y_9 = phi_2_current(1,2) + [q5a(2);q6a(2);q7a(2);q8a(2)];
                    new_z_9 = phi_2_current(1,3) + [q5a(3);q6a(3);q7a(3);q8a(3)];
    
                    new_x_10 = phi_2_current(1,1) + [q4a(1);q3a(1);q7a(1);q8a(1)];
                    new_y_10 = phi_2_current(1,2) + [q4a(2);q3a(2);q7a(2);q8a(2)];
                    new_z_10 = phi_2_current(1,3) + [q4a(3);q3a(3);q7a(3);q8a(3)];
    
                    new_x_11 = phi_2_current(1,1) + [q5a(1);q8a(1);q4a(1);q1a(1)];
                    new_y_11 = phi_2_current(1,2) + [q5a(2);q8a(2);q4a(2);q1a(2)];
                    new_z_11 = phi_2_current(1,3) + [q5a(3);q8a(3);q4a(3);q1a(3)];
    
                    new_x_12 = phi_2_current(1,1) + [q2a(1);q3a(1);q7a(1);q6a(1)];
                    new_y_12 = phi_2_current(1,2) + [q2a(2);q3a(2);q7a(2);q6a(2)];
                    new_z_12 = phi_2_current(1,3) + [q2a(3);q3a(3);q7a(3);q6a(3)];

                    new_x_13 = phi_3_current(1,1) + [r1a(1);r2a(1);r3a(1);r4a(1)];
                    new_y_13 = phi_3_current(1,2) + [r1a(2);r2a(2);r3a(2);r4a(2)];
                    new_z_13 = phi_3_current(1,3) + [r1a(3);r2a(3);r3a(3);r4a(3)];
    
                    new_x_14 = phi_3_current(1,1) + [r1a(1);r2a(1);r6a(1);r5a(1)];
                    new_y_14 = phi_3_current(1,2) + [r1a(2);r2a(2);r6a(2);r5a(2)];
                    new_z_14 = phi_3_current(1,3) + [r1a(3);r2a(3);r6a(3);r5a(3)];
    
                    new_x_15 = phi_3_current(1,1) + [r5a(1);r6a(1);r7a(1);r8a(1)];
                    new_y_15 = phi_3_current(1,2) + [r5a(2);r6a(2);r7a(2);r8a(2)];
                    new_z_15 = phi_3_current(1,3) + [r5a(3);r6a(3);r7a(3);r8a(3)];
    
                    new_x_16 = phi_3_current(1,1) + [r4a(1);r3a(1);r7a(1);r8a(1)];
                    new_y_16 = phi_3_current(1,2) + [r4a(2);r3a(2);r7a(2);r8a(2)];
                    new_z_16 = phi_3_current(1,3) + [r4a(3);r3a(3);r7a(3);r8a(3)];
    
                    new_x_17 = phi_3_current(1,1) + [r5a(1);r8a(1);r4a(1);r1a(1)];
                    new_y_17 = phi_3_current(1,2) + [r5a(2);r8a(2);r4a(2);r1a(2)];
                    new_z_17 = phi_3_current(1,3) + [r5a(3);r8a(3);r4a(3);r1a(3)];
    
                    new_x_18 = phi_3_current(1,1) + [r2a(1);r3a(1);r7a(1);r6a(1)];
                    new_y_18 = phi_3_current(1,2) + [r2a(2);r3a(2);r7a(2);r6a(2)];
                    new_z_18 = phi_3_current(1,3) + [r2a(3);r3a(3);r7a(3);r6a(3)];

                    new_x_19 = phi_4_current(1,1) + [s1a(1);s2a(1);s3a(1);s4a(1)];
                    new_y_19 = phi_4_current(1,2) + [s1a(2);s2a(2);s3a(2);s4a(2)];
                    new_z_19 = phi_4_current(1,3) + [s1a(3);s2a(3);s3a(3);s4a(3)];
    
                    new_x_20 = phi_4_current(1,1) + [s1a(1);s2a(1);s6a(1);s5a(1)];
                    new_y_20 = phi_4_current(1,2) + [s1a(2);s2a(2);s6a(2);s5a(2)];
                    new_z_20 = phi_4_current(1,3) + [s1a(3);s2a(3);s6a(3);s5a(3)];
    
                    new_x_21 = phi_4_current(1,1) + [s5a(1);s6a(1);s7a(1);s8a(1)];
                    new_y_21 = phi_4_current(1,2) + [s5a(2);s6a(2);s7a(2);s8a(2)];
                    new_z_21 = phi_4_current(1,3) + [s5a(3);s6a(3);s7a(3);s8a(3)];
    
                    new_x_22 = phi_4_current(1,1) + [s4a(1);s3a(1);s7a(1);s8a(1)];
                    new_y_22 = phi_4_current(1,2) + [s4a(2);s3a(2);s7a(2);s8a(2)];
                    new_z_22 = phi_4_current(1,3) + [s4a(3);s3a(3);s7a(3);s8a(3)];
    
                    new_x_23 = phi_4_current(1,1) + [s5a(1);s8a(1);s4a(1);s1a(1)];
                    new_y_23 = phi_4_current(1,2) + [s5a(2);s8a(2);s4a(2);s1a(2)];
                    new_z_23 = phi_4_current(1,3) + [s5a(3);s8a(3);s4a(3);s1a(3)];
    
                    new_x_24 = phi_4_current(1,1) + [s2a(1);s3a(1);s7a(1);s6a(1)];
                    new_y_24 = phi_4_current(1,2) + [s2a(2);s3a(2);s7a(2);s6a(2)];
                    new_z_24 = phi_4_current(1,3) + [s2a(3);s3a(3);s7a(3);s6a(3)];

                    plot3(x, y, z, 'k', LineWidth=1.2);
                    hold on
                    plot3(x2, y2, z2, 'k', LineWidth=1.2);
                    plot3(x3, y3, z3, 'k', LineWidth=1.2);
                    plot3(x4, y4, z4, 'k', LineWidth=1.2);
                    plot3(x5, y5, z5, 'k', LineWidth=1.2);
                    plot3(x6, y6, z6, 'k', LineWidth=1.2);
                    plot3(x7, y7, z7, 'k', LineWidth=1.2);
                    plot3(x8, y8, z8, 'k', LineWidth=1.2);

                    fill3(new_x_1,new_y_1, new_z_1,"white")
                    fill3(new_x_2,new_y_2, new_z_2,"white")
                    fill3(new_x_3,new_y_3, new_z_3,"white")
                    fill3(new_x_4,new_y_4, new_z_4,"white")
                    fill3(new_x_5,new_y_5, new_z_5,"white")
                    fill3(new_x_6,new_y_6, new_z_6,"white")
                    fill3(new_x_7,new_y_7, new_z_7,"white")
                    fill3(new_x_8,new_y_8, new_z_8,"white")
                    fill3(new_x_9,new_y_9, new_z_9,"white")
                    fill3(new_x_10,new_y_10, new_z_10,"white")
                    fill3(new_x_11,new_y_11, new_z_11,"white")
                    fill3(new_x_12,new_y_12, new_z_12,"white")
                    fill3(new_x_13,new_y_13, new_z_13,"white")
                    fill3(new_x_14,new_y_14, new_z_14,"white")
                    fill3(new_x_15,new_y_15, new_z_15,"white")
                    fill3(new_x_16,new_y_16, new_z_16,"white")
                    fill3(new_x_17,new_y_17, new_z_17,"white")
                    fill3(new_x_18,new_y_18, new_z_18,"white")
                    fill3(new_x_19,new_y_19, new_z_19,"white")
                    fill3(new_x_20,new_y_20, new_z_20,"white")
                    fill3(new_x_21,new_y_21, new_z_21,"white")
                    fill3(new_x_22,new_y_22, new_z_22,"white")
                    fill3(new_x_23,new_y_23, new_z_23,"white")
                    fill3(new_x_24,new_y_24, new_z_24,"white")

                    plot3(phi_1_current(1,1) + x12(1), phi_1_current(1,2) + x12(2), phi_1_current(1,3)+ x12(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);
                    plot3(phi_1_current(1,1) + x14(1), phi_1_current(1,2) + x14(2), phi_1_current(1,3)+ x14(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);
                    plot3(phi_3_current(1,1) + x32(1), phi_3_current(1,2) + x32(2), phi_3_current(1,3)+ x32(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);
                    plot3(phi_3_current(1,1) + x34(1), phi_3_current(1,2) + x34(2), phi_3_current(1,3)+ x34(3), 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', "black", 'MarkerFaceColor', grayColor);


                    % fill3(x, y, z, 'r')
                    % fill3(x2, y2, z2, 'r')
                    % fill3(x3, y3, z3, 'r')
                    % fill3(x4, y4, z4, 'r')
                    % fill3(x5, y5, z5, 'r')
                    % fill3(x6, y6, z6, 'r')
                    % fill3(x7, y7, z7, 'r')
                    % fill3(x8, y8, z8, 'r')
                    view(136, 23)

                    drawnow

                    %print(gcf, strcat('snapshot_', num2str(j)), '-dpng')

                end

        end

    end

end
