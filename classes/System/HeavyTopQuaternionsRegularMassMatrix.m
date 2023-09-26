%% Class: RigidBodyRotatingQuaternions
%
% The rotation of a rigid body in terms of unit-quaternions.
%

classdef HeavyTopQuaternionsRegularMassMatrix < System

    properties

    end

    methods

        function self = HeavyTopQuaternionsRegularMassMatrix(CONFIG)

            self.mCONSTRAINTS = 1;
            self.nBODIES = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.nDOF = 4;
            self.MASS_MAT = [];
            self.EXT_ACC = CONFIG.EXT_ACC;

           % parameters
            rho = 2700; % mass density
            H = 0.1; % length of the gyro top
            R = 0.05; % radius of the gyro top
            L = 3 * H / 4; % location of center of mass along sym. axis
            MASS = rho * pi * R^2 * H / 3; % total mass of the gyro top
            J1 = 3 / 80 * MASS * (4 * R^2 + H^2)+MASS*L^2; % inertia moment w.r.t. d1-axis (J1 = J2)
            J2 = J1;
            J3 = 3 / 10 * MASS * R^2; % inertia moment w.r.t. d3-axis (sym. axis)
            self.GEOM = [J1,J2,J3,L]; % principal inertia components
            
            self.DISS_MAT = zeros(4,4);
            self.nPotentialInvariants = 0;
            self.nKineticInvariants = 1;
            self.nConstraintInvariants = 1;
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
            M = 4*G_q'*inertia_tensor*G_q;

            % regular mass matrix
            M4 = M + 2*trace(inertia_tensor)*(q*q');

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
        
        % gradient of kinetic energy w.r.t. position formulated in momentum
        function Dq_T = kinetic_energy_gradient_from_momentum(self, q, p)
            
            %extract vector and scalar part form quaternion
            p_vec = p(2:4);
            p_scalar = p(1);

            %skew-sym matrix corresponding to vector part
            p_hat = [0, -p_vec(3), p_vec(2);
                    p_vec(3), 0, -p_vec(1);
                    -p_vec(2), p_vec(1), 0];
            
            % transformation matrix
            G_p = [-p_vec, p_scalar*eye(3) - p_hat];

            Ql_p = [p, G_p'];

            % classical inertia tensor
            inertia_tensor = diag(self.GEOM(1:3));
            extended_inertia_tensor = [1/2*trace(inertia_tensor), zeros(3,1)';
                                       zeros(1,3)', inertia_tensor];
            inv_ext_inertia_tensor = eye(size(extended_inertia_tensor)) / extended_inertia_tensor;

            Dq_T = 1/4 * Ql_p * inv_ext_inertia_tensor * Ql_p' * q;

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
            L = self.GEOM(4);
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
            L = self.GEOM(4);
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
            V_ext = -self.MASS * self.EXT_ACC(3) * L * [0;0;1]' * R_q * [0;0;1];

        end

        function DV_ext = external_potential_gradient(self, q)

            % External potential
            L = self.GEOM(4);
            %extract vector and scalar part form quaternion
            q_vec = q(2:4);
            q_scalar = q(1);

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

            DV_ext = 2*self.MASS*self.EXT_ACC(3)*L*e_bar*E_q'*e_3;

        end

        function D2V_ext = external_potential_hessian(~, ~)

            D2V_ext = zeros(4,4);

        end

        function V_int = internal_potential(~, ~)

            V_int = 0;

        end

        function DV_int = internal_potential_gradient(~, ~)

            DV_int = [0;0; 0; 0];

        end

        function D2V_int = internal_potential_hessian(~, ~)

            D2V_int = diag([0; 0; 0; 0]);

        end


        %% Constraint on position level

        function g = constraint(~, q)

            if size(q,1) == 1 && size(q,2) == 4
                q = q';
            end
           
            g=1/2*(q'*q -1);

        end

        function Dg = constraint_gradient(~, q)

            if size(q,1) == 1 && size(q,2) == 4
                q = q';
            end
           
            Dg=q';

        end

        function D2g = constraint_hessian(~, ~, ~)

             D2g=eye(4);

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
        function omega = kinetic_energy_invariant(~, q, p, ~)
            
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
            
            Ql_q = [q, G_q'];

            omega = Ql_q'*p;

        end

            % gradient of potential invariants w.r.t. q
        function DomegaDq = kinetic_energy_invariant_gradient_q(~, ~, p, ~)
            %extract vector and scalar part form quaternion
            p_vec = p(2:4);
            p_scalar = p(1);

            %skew-sym matrix corresponding to vector part
            p_hat = [0, -p_vec(3), p_vec(2);
                    p_vec(3), 0, -p_vec(1);
                    -p_vec(2), p_vec(1), 0];
            
            % transformation matrix
            G_p = [-p_vec, p_scalar*eye(3) - p_hat];

            DomegaDq = [p'; 
                       -G_p]; %=Q_r(p)*I_quer;

        end

         % gradient of kinetic energy invariants w.r.t. p
        function DomegaDp = kinetic_energy_invariant_gradient_p(~, q, ~, ~)
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

            DomegaDp = [q';
                        G_q];
        end

        % kinetic energy computed with the invariant
        function Ts = kinetic_energy_from_invariant(self, omega, ~)
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM(1:3));
            extended_inertia_tensor = [1/2*trace(inertia_tensor), zeros(3,1)';
                                       zeros(1,3)', inertia_tensor];
            inv_ext_inertia_tensor = eye(size(extended_inertia_tensor)) / extended_inertia_tensor;
            Ts = 1/8 *omega'*inv_ext_inertia_tensor*omega;
        end

        % gradient of kinetic energy w.r.t. the invariant
        function DTsDpi = kinetic_energy_gradient_from_invariant(self, omega, ~)
          
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM(1:3));
            extended_inertia_tensor = [1/2*trace(inertia_tensor), zeros(3,1)';
                                       zeros(1,3)', inertia_tensor];
            inv_ext_inertia_tensor = eye(size(extended_inertia_tensor)) / extended_inertia_tensor;

            DTsDpi = 1/4*inv_ext_inertia_tensor*omega;
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
       
       if size(q,1) == 1 && size(q,2) == 4
            q = q';
       end

       zeta = q'*q;

    end

    % gradient of the invariant of the position constraint w.r.t. q
    function DzetaDq = constraint_invariant_gradient(~, q, ~)
        if size(q,1) == 1 && size(q,2) == 4
            q = q';
        end
        DzetaDq = 2*q';
    end

    % gradient of the invariant of the position constraint w.r.t. q
    function D2zetaDq2 = constraint_invariant_hessian(~, ~, ~)

      D2zetaDq2 = 2*eye(4);
    end

    % position constrained computed with its invariant
      function gs = constraint_from_invariant(~, zeta, ~)

           gs = 1/2*(zeta -1);
      end

    % gradient of position constrained w.r.t. its invariant
    function Dgs = constraint_gradient_from_invariant(~, ~, ~)

        Dgs = 1/2;

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
                q1 = this_simulation.z(:, 1:DIM);
                q2 = this_simulation.z(:, DIM+1:2*DIM);
                q3 = this_simulation.z(:, 2*DIM+1:3*DIM);
                q4 = this_simulation.z(:, 3*DIM+1:4*DIM);
                xmin = min([min(q1(:, 1)), min(q2(:, 1)), min(q3(:, 1)), min(q4(:, 1))]);
                ymin = min([min(q1(:, 2)), min(q2(:, 2)), min(q3(:, 2)), min(q4(:, 2))]);
                zmin = min([min(q1(:, 3)), min(q2(:, 3)), min(q3(:, 3)), min(q4(:, 3))]);
                xmax = max([max(q1(:, 1)), max(q2(:, 1)), max(q3(:, 1)), max(q4(:, 1))]);
                ymax = max([max(q1(:, 2)), max(q2(:, 2)), max(q3(:, 2)), max(q4(:, 2))]);
                zmax = max([max(q1(:, 3)), max(q2(:, 3)), max(q3(:, 3)), max(q4(:, 3))]);
                NT = size(q1, 1);

                axis equal
                %axis([xmin, xmax, ymin, ymax, zmin, zmax]);
                %xlabel('x');
                %ylabel('y');
                %zlabel('z');
                %grid on;
                axis off

                xa1 = q1(1, 1);
                xa2 = q2(1, 1);
                xa3 = q3(1, 1);
                xa4 = q4(1, 1);
                ya1 = q1(1, 2);
                ya2 = q2(1, 2);
                ya3 = q3(1, 2);
                ya4 = q4(1, 2);

                if DIM == 3
                    za1 = q1(1, 3);
                    za2 = q2(1, 3);
                    za3 = q3(1, 3);
                    za4 = q4(1, 3);
                else
                    za1 = 0;
                    za2 = 0;
                    za3 = 0;
                    za4 = 0;
                end

                for j = 1:NT

                    cla(fig);
                    hold on

                    %% Current position
                    x1 = q1(j, 1);
                    x2 = q2(j, 1);
                    x3 = q3(j, 1);
                    x4 = q4(j, 1);

                    y1 = q1(j, 2);
                    y2 = q2(j, 2);
                    y3 = q3(j, 2);
                    y4 = q4(j, 2);

                    if DIM == 3
                        z1 = q1(j, 3);
                        z2 = q2(j, 3);
                        z3 = q3(j, 3);
                        z4 = q4(j, 3);
                    else
                        z1 = 0;
                        z2 = 0;
                        z3 = 0;
                        z4 = 0;
                    end

                    %% Reference sphere
                    %                     plot3(xa1, ya1, za1, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    %                     hold on
                    %                     plot3(xa2, ya2, za2, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    %                     hold on
                    %                     plot3(xa3, ya3, za3, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    %                     hold on
                    %                     plot3(xa4, ya4, za4, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    %                     hold on

                    %% Reference constraint
                    %                     xx3 = [xa1; xa2];
                    %                     yy3 = [ya1; ya2];
                    %                     zz3 = [za1; za2];
                    %                     xxx3 = [xa3; xa4];
                    %                     yyy3 = [ya3; ya4];
                    %                     zzz3 = [za3; za4];
                    %                     plot3(xx3, yy3, zz3, 'k', 'LineWidth', 2);
                    %                     hold on
                    %                     plot3(xxx3, yyy3, zzz3, 'k', 'LineWidth', 2);

                    %% Reference position of the springs
                    %                     xx_3 = [xa1; xa3];
                    %                     yy_3 = [ya1; ya3];
                    %                     zz_3 = [za1; za3];
                    %                     xxxx_3 = [xa2; xa4];
                    %                     yyyy_3 = [ya2; ya4];
                    %                     zzzz_3 = [za2; za4];
                    %
                    %                     tmp  = gfx_springelement([xa1,ya1,za1],[xa3,ya3,za3],1,0.1,5);
                    %                     tmp2 = gfx_springelement([xa2,ya2,za2],[xa4,ya4,za4],1,0.1,5);
                    %
                    %                     %plot3(xx_3, yy_3, zz_3, 'k--');
                    %                     plot3(tmp(1,:),tmp(2,:),tmp(3,:),'LineWidth',1.5);
                    %                     hold on
                    %                     %plot3(xxxx_3, yyyy_3, zzzz_3, 'k--');
                    %                     plot3(tmp2(1,:),tmp2(2,:),tmp2(3,:),'LineWidth',1.5);

                    %% current position of the mass
                    %hold on
                    %if DIM == 3
                    %    plot3(q1(1:j, 1), q1(1:j, 2), q1(1:j, 3), 'k');
                    %    plot3(q2(1:j, 1), q2(1:j, 2), q2(1:j, 3), 'k');
                    %    plot3(q3(1:j, 1), q3(1:j, 2), q3(1:j, 3), 'k');
                    %    plot3(q4(1:j, 1), q4(1:j, 2), q4(1:j, 3), 'k');
                    %else
                    %    plot3(q1(1:j, 1), q1(1:j, 2), zeros(j, 1), 'k');
                    %    plot3(q2(1:j, 1), q2(1:j, 2), zeros(j, 1), 'k');
                    %    plot3(q3(1:j, 1), q3(1:j, 2), zeros(j, 1), 'k');
                    %    plot3(q4(1:j, 1), q4(1:j, 2), zeros(j, 1), 'k');
                    %end
                    plot3(x1, y1, z1, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
                    plot3(x2, y2, z2, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
                    plot3(x3, y3, z3, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
                    plot3(x4, y4, z4, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');

                    %grid on

                    %% current position of the constraint
                    x_3 = [x1; x2];
                    y_3 = [y1; y2];
                    z_3 = [z1; z2];
                    xxx3 = [x3; x4];
                    yyy3 = [y3; y4];
                    zzz3 = [z3; z4];
                    plot3(x_3, y_3, z_3, 'k', 'linewidth', 2);
                    plot3(xxx3, yyy3, zzz3, 'k', 'linewidth', 2);

                    %% Current position of the springs
                    xx3 = [x1; x3];
                    yy3 = [y1; y3];
                    zz3 = [z1; z3];
                    xxxx3 = [x2; x4];
                    yyyy3 = [y2; y4];
                    zzzz3 = [z2; z4];

                    tmp3 = gfx_springelement([x1, y1, z1], [x3, y3, z3], 1, 0.05, 7);
                    tmp4 = gfx_springelement([x2, y2, z2], [x4, y4, z4], 1, 0.05, 7);

                    %plot3(xx3, yy3, zz3, 'k--');
                    plot3(tmp3(1, :), tmp3(2, :), tmp3(3, :), 'k', 'LineWidth', 1.5);
                    hold on
                    plot3(tmp4(1, :), tmp4(2, :), tmp4(3, :), 'k', 'LineWidth', 1.5);
                    %plot3(xxxx3, yyyy3, zzzz3, 'k--');


                    if DIM == 2
                        view(0, 90)
                    else
                        view(136, 23)
                    end

                    drawnow

                    print(gcf, strcat('snapshot_', num2str(j)), '-dpng')

                end

        end

    end

end
