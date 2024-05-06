%% Class: PendulumQuaternionsRot
%
% The rotation of a rigid body in terms of unit-quaternions.
%

classdef PendulumQuaternionsRigidBody < System

    properties

    end

    methods

        function self = PendulumQuaternionsRigidBody(CONFIG)

            self.mCONSTRAINTS = 2;
            self.nBODIES = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.nDOF = 4;
            self.MASS_MAT = [];
            self.EXT_ACC = CONFIG.EXT_ACC;

           % parameters
            self.GEOM = [self.MASS*5^2,self.MASS*5^2,0,5]; % mass / length
            
            self.DISS_MAT = zeros(4,4);
            self.nPotentialInvariants = 0;
            self.nKineticInvariants = 0;
            self.nConstraintInvariants = 2;
            self.mMixedQuantities = 0;
            self.isCyclicCoordinate = [false;false;false;false];

        end
        
        function M4 = get_mass_matrix(self, q)
            
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
            
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM(1:3));
            
            % singular mass matrix
            M4 = 4*G_q'*inertia_tensor*G_q;

        end

        function Dq_T = kinetic_energy_gradient_from_velocity(self, q, v)
            
            %extract vector and scalar part form quaternion
            v_vec = v(2:4);
            v_scalar = v(1);

            %skew-sym matrix corresponding to vector part
            v_hat = [0, -v_vec(3), v_vec(2);
                    v_vec(3), 0, -v_vec(1);
                    -v_vec(2), v_vec(1), 0];
            
            % transformation matrix
            G_v = [-v_vec, v_scalar*eye(3) - v_hat];
            
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM(1:3));
            
            % singular inertia matrix in v
            M_4_hat = 4*G_v'*inertia_tensor*G_v;
            
            % partial derivativa of kinetic energy w.r.t. q(uat)
            Dq_T = M_4_hat * q;

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(~, ~, ~)
            
            Dq_T = [];

        end


        function L = get_cartesian_angular_momentum_from_momentum(~, q, p)            
            
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
            
            L = 1/2 * E_q * p ;

        end

        function x = get_cartesian_coordinates_center_of_mass(self,q)
            
            % External potential
            L = self.GEOM(end);

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
            E_q = [-q_vec, q_scalar*eye(3) + q_hat];
            R_q = E_q * G_q';
            d3 = R_q*[0;0;1];
            x = L * d3;

        end

        %% Potential functions

        function V_ext = external_potential(self, q)
            % External potential
            m = self.MASS;
            l = self.GEOM(end);
            g = self.EXT_ACC;

            q_vec = q(2:4);
            q_scalar = q(1);

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
            L = l;
            % transformation matrix
            G_q = [-q_vec, q_scalar*eye(3) - q_hat];
            E_q = [-q_vec, q_scalar*eye(3) + q_hat];
            R_q = E_q * G_q';
            V_ext = -m * g * L * [0;0;1]' * R_q * [0;0;1];


        end

        function DV_ext = external_potential_gradient(self, q)
            
            %extract vector and scalar part form quaternion
            q_scalar = q(1);
            q_vec = q(2:4);
            m = self.MASS;
            l = self.GEOM(end);
            g = self.EXT_ACC;

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
            % transformation matrix
            E_q = [-q_vec, q_scalar*eye(3) + q_hat];
            
            e_3 = [0;0;1];
            e_hat = [0, -e_3(3), e_3(2);
                    e_3(3), 0, -e_3(1);
                    -e_3(2), e_3(1), 0];
            e_bar = [0, -e_3';
                     e_3, -e_hat];
            L = l;

            DV_ext = 2*m*g*L*e_bar*E_q'*e_3;

        end

        function D2V_ext = external_potential_hessian(~, ~)

            D2V_ext = zeros(4,4);

        end

        function V_int = internal_potential(~, ~)

            V_int = 0;

        end

        function DV_int = internal_potential_gradient(~, ~)

            DV_int = [0; 0; 0; 0];

        end

        function D2V_int = internal_potential_hessian(~, ~)

            D2V_int = diag([0; 0; 0; 0]);

        end


        %% Constraint on position level

        function g = constraint(~, q)
           
            g1=1/2*(q'*q - 1);
            %g2= q(1)*q(4)-q(2)*q(3);
            g2 = q(1)*q(2)+q(3)*q(4);
            g = [g1;g2];
            %g = g1;
        end

        function Dg = constraint_gradient(~, q)
           
            Dg1 = q';
            %Dg2 = [q(4) -q(3) -q(2) q(1)];
            Dg2 = [q(2) q(1) q(4) q(3)];
            Dg = [Dg1; Dg2];
            %Dg = Dg1;

        end

        function D2g = constraint_hessian(~, ~, i)
             if i ==1
                 D2g=eye(4);
             elseif i ==2
                 D2g= [0 0 0 1;
                       0 0 -1 0;
                       0 -1 0 1;
                       1 0 0 0];
             end
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
        function omega = kinetic_energy_invariant(~, q, v, ~)
            
            omega = [];

        end

            % gradient of potential invariants w.r.t. q
        function DomegaDq = kinetic_energy_invariant_gradient_q(~, ~, v, ~)

            % derivative of minimal velocity w.r.t quaternion q
            DomegaDq = [];

        end

         % gradient of kinetic energy invariants w.r.t. v
        function Dv_minDv = kinetic_energy_invariant_gradient_v(~, q, ~, ~)

           Dv_minDv = [];%2*H';
        end

        % internal potential computed with the invariant
        function Ts = kinetic_energy_from_invariant(self, v_min, ~)

            Ts = [];
        end

        % gradient of internal potential w.r.t. the invariant
        function DTsDpi = kinetic_energy_gradient_from_invariant(self, v_min, q)
          
            DTsDpi = [];
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
    function zeta = constraint_invariant(~, q, j)
        if j == 1
           zeta = q'*q;
        elseif j == 2
            zeta = q(1)*q(2)+q(3)*q(4);
        end
    end

    % gradient of the invariant of the position constraint w.r.t. q
    function DzetaDq = constraint_invariant_gradient(~, q, j)
        if j == 1
            DzetaDq = 2*q';
        elseif j ==2
            DzetaDq = [q(2) q(1) q(4) q(3)];
        end
    end

    % gradient of the invariant of the position constraint w.r.t. q
    function D2zetaDq2 = constraint_invariant_hessian(~, ~, j)
        if j ==1
            D2zetaDq2 = 2*eye(4);
        elseif j ==2
             D2zetaDq2 = [0 1 0 0;
                        1 0 0 0;
                        0 0 0 1;
                        0 0 1 0];
        end
    end

    % position constrained computed with its invariant
      function gs = constraint_from_invariant(~, zeta, j)
           if j == 1
               gs = 1/2*(zeta -1);
           elseif j == 2
               gs = zeta;
           end
      end

    % gradient of position constrained w.r.t. its invariant
    function Dgs = constraint_gradient_from_invariant(~, ~, j)
        if j == 1
            Dgs = 1/2;
        elseif j ==2
            Dgs = 1;
        end

    end

     function analyzed_quantity = hconvergence_set(~, this_simulation)

        if strcmp(this_simulation.CONV_QUANTITY,'q')

            q = this_simulation.z(end,1:4)';
            %extract vector and scalar part form quaternion
            q_vec = q(2:4);
            q_scalar = q(1);

            %skew-sym matrix corresponding to vector part
            q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
            % transformation matrix
            G_q = [-q_vec, q_scalar*eye(3) - q_hat];
            E_q = [-q_vec, q_scalar*eye(3) + q_hat];

            R_q = E_q*G_q';

            analyzed_quantity = R_q(:); %rotation matrix at t_end
      

        elseif strcmp(this_simulation.CONV_QUANTITY,'p')
            error('not available.')
        elseif strcmp(this_simulation.CONV_QUANTITY,'lambda')
            error('not available.')
        else
            error('quantity not yet implemented for convergence analysis.')
        end

    end


    function [reference_solution, this_simulation] = hconvergence_reference(~, this_simulation, analyzed_quantity)

        reference_solution = analyzed_quantity(:, end, end); %position
        this_simulation.matrix_error_analysis = true;
      
    end

        %% Animation method
         function give_animation(self, fig, this_simulation)

                DIM = self.DIM;
                quat = this_simulation.z(:, 1:4);
                l = self.GEOM(end);
                NT = size(quat, 1);

                axis equal
                axis([-1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l]);
                xlabel('x');
                ylabel('y');
                zlabel('z');
                grid on;
                
                qa = get_cartesian_coordinates_center_of_mass(self,quat(1,:))';

                xa = qa(1, 1);
                ya = qa(1, 2);
                if DIM == 3
                    za = qa(1, 3);
                else
                    za = 0;
                end

                for j = 1:NT
                    q(j,:) = get_cartesian_coordinates_center_of_mass(self,quat(j,:))';
                end

                for j = 1:NT

                    cla(fig);
                    hold on

                    %% Current position
                    %q = get_cartesian_coordinates_center_of_mass(self,quat(j,:))';
                    x = q(j,1);
                    y = q(j,2);
                    if DIM == 3
                        z = q(j,3);
                    else
                        z = 0;
                    end

                    %% Reference sphere
                    plot3(xa, ya, za, 'mo', 'MarkerSize', 40);%, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
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
