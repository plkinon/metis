%% Class: RigidBodyRotatingQuaternions
%
% The rotation of a rigid body in terms of unit-quaternions.
%

classdef RigidBodyRotatingQuaternions < System

    properties

    end

    methods

        function self = RigidBodyRotatingQuaternions(CONFIG)

            self.mCONSTRAINTS = 1;
            self.nBODIES = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.nDOF = 4;
            self.MASS_MAT = [];
            self.EXT_ACC = CONFIG.EXT_ACC;

           % parameters
            self.GEOM = [6,8,3]; % principal inertia components
            
            self.DISS_MAT = zeros(4,4);
            self.nPotentialInvariants = 0;
            self.nKineticInvariants = 1;
            self.nConstraintInvariants = 1;
            self.mMixedQuantities = 0;
            self.isCyclicCoordinate = [false;false;false;false];

        end
        
        function M = get_mass_matrix(self, q)
            
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
            inertia_tensor = diag(self.GEOM);
            
            % singular mass matrix
            M = 4*G_q'*inertia_tensor*G_q;

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
            inertia_tensor = diag(self.GEOM);
            
            % singular inertia matrix in v
            M_4_hat = 4*G_v'*inertia_tensor*G_v;
            
            % partial derivativa of kinetic energy w.r.t. q(uat)
            Dq_T = M_4_hat * q;

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(~, ~, ~)
            
            Dq_T = [0; 0; 0];

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

        %% Potential functions

        function V_ext = external_potential(~, ~)
            % External potential
            V_ext = 0;

        end

        function DV_ext = external_potential_gradient(~, ~)

            DV_ext = [0; 0; 0; 0];

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
           
            g=1/2*(q'*q -1);

        end

        function Dg = constraint_gradient(~, q)
           
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
        function omega = kinetic_energy_invariant(~, q, v, ~)
            
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

            omega = 2*G_q*v;

        end

            % gradient of potential invariants w.r.t. q
        function DomegaDq = kinetic_energy_invariant_gradient_q(~, ~, v, ~)
            %extract vector and scalar part form quaternion
            v_vec = v(2:4);
            v_scalar = v(1);

            %skew-sym matrix corresponding to vector part
            v_hat = [0, -v_vec(3), v_vec(2);
                    v_vec(3), 0, -v_vec(1);
                    -v_vec(2), v_vec(1), 0];
            
            % transformation matrix
            G_v = [-v_vec, v_scalar*eye(3) - v_hat];
            DomegaDq = -2*G_v;
        end

         % gradient of kinetic energy invariants w.r.t. v
        function DomegaDv = kinetic_energy_invariant_gradient_v(~, q, ~, ~)
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
            DomegaDv = 2*G_q;
        end

        % internal potential computed with the invariant
        function Ts = kinetic_energy_from_invariant(self, omega, ~)
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM);

            Ts = 1/2 *omega'*inertia_tensor*omega;
        end

        % gradient of internal potential w.r.t. the invariant
        function DTsDpi = kinetic_energy_gradient_from_invariant(self, omega, ~)
          
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM);

            DTsDpi = inertia_tensor*omega;
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

    end

    % gradient of the invariant of the position constraint w.r.t. q
    function DzetaDq = constraint_invariant_gradient(~, q, ~)

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

    function qn1 = get_coordinate_update(~,theta,qn)

        abs_theta = norm(theta);
        if abs_theta == 0
            expS3 = [1;0;0;0];
        else
            expS3 = [cos(1/2*abs_theta); 1/abs_theta * sin(1/2*abs_theta) * theta];
        end
        

        exp_vec = expS3(2:4);
        exp_scalar = expS3(1);

        %skew-sym matrix corresponding to vector part
        exp_hat = [0, -exp_vec(3), exp_vec(2);
                exp_vec(3), 0, -exp_vec(1);
                -exp_vec(2), exp_vec(1), 0];
        
        % transformation matrix
        G_exp = [-exp_vec, exp_scalar*eye(3) - exp_hat];

        Ql_exp = [expS3, G_exp'];

        qn1 = Ql_exp*qn;

    end

    function P = get_null_space_matrix(~,q)

        q_vec = q(2:4);
        q_scalar = q(1);

        %skew-sym matrix corresponding to vector part
        q_hat = [0, -q_vec(3), q_vec(2);
                q_vec(3), 0, -q_vec(1);
                -q_vec(2), q_vec(1), 0];
        
        % transformation matrix
        P = [-q_vec, q_scalar*eye(3) - q_hat];

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
            function give_animation(~, fig, this_simulation)

                q = this_simulation.z(:, 1:4);
                qa = q(1,:);
                NT = size(q, 1);

                %extract vector and scalar part form quaternion
                q_vec = qa(2:4)';
                q_scalar = qa(1);

                q_hat = [0, -q_vec(3), q_vec(2);
                    q_vec(3), 0, -q_vec(1);
                    -q_vec(2), q_vec(1), 0];
            
                % transformation matrix
                G_q = [-q_vec, q_scalar*eye(3) - q_hat];
                E_q = [-q_vec, q_scalar*eye(3) + q_hat];
                R_q = E_q*G_q';

                b=6;
                h=8;
                t=3;

                P1 = [t/2,b/2,-h/2]';
                P2 = [-t/2,b/2,-h/2]';
                P3 = [-t/2,-b/2,-h/2]';
                P4 = [t/2,-b/2,-h/2]';
                P5 = [t/2,b/2,h/2]';
                P6 = [-t/2,b/2,h/2]';
                P7 = [-t/2,-b/2,h/2]';
                P8 = [t/2,-b/2,h/2]';

                p1a = R_q*P1;
                p2a = R_q*P2;
                p3a = R_q*P3;
                p4a = R_q*P4;
                p5a = R_q*P5;
                p6a = R_q*P6;
                p7a = R_q*P7;
                p8a = R_q*P8;

                xa = [p1a(1);p2a(1);p3a(1);p4a(1);p1a(1);p5a(1);p6a(1);p7a(1);p8a(1);p5a(1)];
                ya = [p1a(2);p2a(2);p3a(2);p4a(2);p1a(2);p5a(2);p6a(2);p7a(2);p8a(2);p5a(2)];
                za = [p1a(3);p2a(3);p3a(3);p4a(3);p1a(3);p5a(3);p6a(3);p7a(3);p8a(3);p5a(3)];
                
                x2a = [p4a(1);p8a(1);p7a(1);p3a(1);p2a(1);p6a(1)];
                y2a = [p4a(2);p8a(2);p7a(2);p3a(2);p2a(2);p6a(2)];
                z2a = [p4a(3);p8a(3);p7a(3);p3a(3);p2a(3);p6a(3)];

                d1 = R_q*[1;0;0];
                d2 = R_q*[0;1;0];
                d3 = R_q*[0;0;1];
                
                % x = [p1a'; p2a'; p3a'; p4a'];
                % y = [p4a'; p3a'; p7a'; p8a'];
                % z = [p2a'; p3a'; p7a'; p5a'];
                
                plot3(xa, ya, za, 'k', LineWidth=1.2);
                hold on
                plot3(x2a, y2a, z2a, 'k', LineWidth=1.2);
                h1 = mArrow3([0 0 0],2*d1', 'facealpha', 0.5, 'color', 'blue', 'stemWidth', 0.04);
                h2 = mArrow3([0 0 0],2*d2', 'facealpha', 0.5, 'color', 'blue', 'stemWidth', 0.04);
                h3 = mArrow3([0 0 0],2*d3', 'facealpha', 0.5, 'color', 'blue', 'stemWidth', 0.04);

                axis equal
                xmin = -5.2;
                xmax = 5.2;
                ymin = -5.2;
                ymax = 5.2;
                zmin =-5.2;
                zmax =5.2;
                axis([xmin, xmax, ymin, ymax, zmin, zmax]);

                grid on;
                axis off

                for j = 1:NT

                    cla(fig);
                    hold on

                    %% Current position
                    q_current = q(j,:);
                    
                    q_vec = q_current(2:4)';
                    q_scalar = q_current(1);

                    q_hat = [0, -q_vec(3), q_vec(2);
                            q_vec(3), 0, -q_vec(1);
                           -q_vec(2), q_vec(1), 0];
            
                    % transformation matrix
                    G_q = [-q_vec, q_scalar*eye(3) - q_hat];
                    E_q = [-q_vec, q_scalar*eye(3) + q_hat];
                    R_q = E_q*G_q';
    
                    p1a = R_q*P1;
                    p2a = R_q*P2;
                    p3a = R_q*P3;
                    p4a = R_q*P4;
                    p5a = R_q*P5;
                    p6a = R_q*P6;
                    p7a = R_q*P7;
                    p8a = R_q*P8;
    
                    x = [p1a(1);p2a(1);p3a(1);p4a(1);p1a(1);p5a(1);p6a(1);p7a(1);p8a(1);p5a(1)];
                    y = [p1a(2);p2a(2);p3a(2);p4a(2);p1a(2);p5a(2);p6a(2);p7a(2);p8a(2);p5a(2)];
                    z = [p1a(3);p2a(3);p3a(3);p4a(3);p1a(3);p5a(3);p6a(3);p7a(3);p8a(3);p5a(3)];
                    
                    x2 = [p4a(1);p8a(1);p7a(1);p3a(1);p2a(1);p6a(1)];
                    y2 = [p4a(2);p8a(2);p7a(2);p3a(2);p2a(2);p6a(2)];
                    z2 = [p4a(3);p8a(3);p7a(3);p3a(3);p2a(3);p6a(3)];
                    
                    d1t = R_q*[1;0;0];
                    d2t = R_q*[0;1;0];
                    d3t = R_q*[0;0;1];
                
                    plot3(x, y, z, 'k', LineWidth=1.2);
                    hold on
                    plot3(x2, y2, z2, 'k', LineWidth=1.2);
                    h1t = mArrow3([0 0 0],2*d1t', 'facealpha', 1.0, 'color','#CC79A7', 'stemWidth', 0.04);
                    h2t = mArrow3([0 0 0],2*d2t', 'facealpha', 1.0, 'color','#CC79A7', 'stemWidth', 0.04);
                    h3t = mArrow3([0 0 0],2*d3t', 'facealpha', 1.0, 'color', '#CC79A7', 'stemWidth', 0.04);


                    view(136, 23)

                    drawnow

                    print(gcf, strcat('snapshot_', num2str(j)), '-dpng')

                end

        end

    end

end
