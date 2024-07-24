%% Class: RigidBodyRotatingQuaternions
%
% The rotation of a rigid body in terms of unit-quaternions.
%

classdef HeavyTopQuaternions < System

    properties

    end

    methods

        function self = HeavyTopQuaternions(CONFIG)

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
            inertia_tensor = diag(self.GEOM(1:3));

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
            inertia_tensor = diag(self.GEOM(1:3));
            
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

        function x = get_cartesian_coordinates_center_of_mass(self,q)
            
            % External potential
            L = self.GEOM(4);

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
            inertia_tensor = diag(self.GEOM(1:3));

            Ts = 1/2 *omega'*inertia_tensor*omega;
        end

        % gradient of internal potential w.r.t. the invariant
        function DTsDpi = kinetic_energy_gradient_from_invariant(self, omega, ~)
          
            % classical inertia tensor
            inertia_tensor = diag(self.GEOM(1:3));

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

     function analyzed_quantity = hconvergence_set(self, this_simulation)

        if strcmp(this_simulation.CONV_QUANTITY,'q')

            q = this_simulation.z(end,1:4)';
            phi = self.get_cartesian_coordinates_center_of_mass(q);

            analyzed_quantity = phi; 
      

        elseif strcmp(this_simulation.CONV_QUANTITY,'p')
            error('not available.')
        elseif strcmp(this_simulation.CONV_QUANTITY,'lambda')
            error('not available.')
        else
            error('quantity not yet implemented for convergence analysis.')
        end

    end


    function [reference_solution, this_simulation] = hconvergence_reference(self, this_simulation, ~)

        L = self.GEOM(4);
        theta_0 = pi/3;
        omega_p = 10;
        t=this_simulation.t(end);
        reference_solution = [L*sin(theta_0)*sin(omega_p*t);
                              -L*sin(theta_0)*cos(omega_p*t);
                              L*cos(theta_0)];%position
        this_simulation.matrix_error_analysis = false;
      
    end

        %% Animation method
            function give_animation(self, fig, this_simulation)

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
                
                ampl_factor = 10;
                a=ampl_factor*0.1;
                r = a/2;
                l=3/4*a;
                
                P1 = [0,0,a]';
                P2 = [0,0,l]';

                p1a = R_q*P1;
                p2a = R_q*P2;
                
                d1 = R_q*[0.8;0;0];
                d2 = R_q*[0;0.8;0];
                d3 = R_q*[0;0;0.8];

                h1 = mArrow3(p2a,p2a'+d1', 'facealpha', 0.5, 'color', [70/255 100/255 170/255], 'stemWidth', 0.015);
                h2 = mArrow3(p2a,p2a'+d2', 'facealpha', 0.5, 'color', [70/255 100/255 170/255], 'stemWidth', 0.015);
                h3 = mArrow3(p2a,p2a'+d3', 'facealpha', 0.5, 'color', [70/255 100/255 170/255], 'stemWidth', 0.015);
                hold on
                % [X,Y,Z]=cylinder([0 1],1000);
                % M=makehgtform('translate',[-1,-1,-1],'xrotate',pi/4,'yrotate',pi/4);
                % h=surf(X,Y,Z,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.3);  
                
                d(1,:) = zeros(3,1); 
                d(2,:) = p1a;
                [X3,Y3,Z3,~,~,~] = cone(atan(r/a)/(2*pi)*360,d,a);
                surf(X3,Y3,Z3,'LineStyle','none','FaceColor',[0 150/255 130/255],'FaceAlpha',0.8);

                axis equal
                xmin = -1.4;
                xmax = 1.4;
                ymin = -1.4;
                ymax = 1.4;
                zmin =-0.3;
                zmax =1.2;
                axis([xmin, xmax, ymin, ymax, zmin, zmax]);

                grid on;
                axis off;

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

                    d1 = R_q*[0.8;0;0];
                    d2 = R_q*[0;0.8;0];
                    d3 = R_q*[0;0;0.8];

                    h1 = mArrow3(p2a,p2a'+d1', 'facealpha', 0.8, 'color', [70/255 100/255 170/255], 'stemWidth', 0.015);
                    h2 = mArrow3(p2a,p2a'+d2', 'facealpha', 0.8, 'color', [70/255 100/255 170/255], 'stemWidth', 0.015);
                    h3 = mArrow3(p2a,p2a'+d3', 'facealpha', 0.8, 'color', [70/255 100/255 170/255], 'stemWidth', 0.015);
                    hold on
                    % [X,Y,Z]=cylinder([0 1],1000);
                    % M=makehgtform('translate',[-1,-1,-1],'xrotate',pi/4,'yrotate',pi/4);
                    % h=surf(X,Y,Z,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.3);  
                    
                    d(1,:) = zeros(3,1); 
                    d(2,:) = p1a;
                    [X3,Y3,Z3,~,~,~] = cone(atan(r/a)/(2*pi)*360,d,a);
                    surf(X3,Y3,Z3,'LineStyle','none','FaceColor',[0 150/255 130/255],'FaceAlpha',0.8, 'EdgeColor',[0 0 0], 'LineWidth',1);

                    view(136, 23)
                    drawnow
                    print(gcf, strcat('snapshot_', num2str(j)), '-dpng')

                end

        end

    end

end
