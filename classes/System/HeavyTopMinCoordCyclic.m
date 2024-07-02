%% Class: PendulumMinCoord
%
% Spherical pendulum in minimal coordinates (two angles).
%

classdef HeavyTopMinCoordCyclic < System

    properties

    end

    methods

        function self = HeavyTopMinCoordCyclic(CONFIG)

            self.mCONSTRAINTS = 0;
            self.nBODIES = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.nDOF = 3;
            self.mMixedQuantities = 0;
            self.MASS_MAT = [];
            self.EXT_ACC = CONFIG.EXT_ACC;

            % Geometric parameters (correspond to symmetric cone)
            a = 0.1; % total height of the cone
            r = a / 2; % top radius of the cone
            l = 3 * a / 4; % location of center of mass along symmetry axis

            % Principle moments of inertia
            J1 = 3 / 80 * self.MASS * (4 * r^2 + a^2) + self.MASS * l^2;
            J2 = J1;
            J3 = 3 / 10 * self.MASS * r^2;

            self.GEOM = [J1, J2, J3, l];
            self.nPotentialInvariants = 0;
            self.mMixedQuantities = 0;
            
            self.DISS_MAT = zeros(3,3);

            self.isCyclicCoordinate = [false;true;true];
        

        end
        
        function M = get_mass_matrix(self, q)
            
            theta = q(1);
            J_1 = self.GEOM(1);
            J_3 = self.GEOM(3);
            M = [J_1 0 0 ;
                0 J_3 J_3*cos(theta);
                0 J_3*cos(theta) J_3*cos(theta)^2+J_1*sin(theta)^2];

        end

        function M = get_mass_matrix_cyclic(self, x)
            
            theta = x(1);
            J_1 = self.GEOM(1);
            J_3 = self.GEOM(3);
            M = [J_1 0 0 ;
                0 J_3 J_3*cos(theta);
                0 J_3*cos(theta) J_3*cos(theta)^2+J_1*sin(theta)^2];

        end

        function r = get_cartesian_coordinates(self, q)

            theta = q(1);
            phi = q(2);
            l = self.GEOM(1);
            r = [l*sin(theta)*cos(phi); l*sin(theta)*sin(phi); l*cos(theta)];
            
        end

        function r_dot = get_cartesian_velocities(self, q, v)
            
            % theta = q(1);
            % phi = q(2);
            % v_theta = v(1);
            % v_phi = v(2);
            % l = self.GEOM(1);
            % r_dot = [l*cos(theta)*cos(phi)*v_theta - l*sin(theta)*sin(phi)*v_phi; l*cos(theta)*sin(phi)*v_theta + l*sin(theta)*cos(phi)*v_phi; l*sin(theta)*v_theta];
            
            error('Not implemented.')

        end

        function r_dot = get_cartesian_velocities_from_momentum(self, q, p)
            
            % theta = q(1);
            % phi = q(2);
            % l = self.GEOM(1);
            % m = self.MASS;
            % p_theta = p(1);
            % v_theta = 1/(m*l^2)*p_theta;
            % p_phi = p(2);
            % v_phi = 1/(m*l^2*sin(theta)^2)*p_phi;
            % r_dot = [l*cos(theta)*cos(phi)*v_theta - l*sin(theta)*sin(phi)*v_phi; 
            %          l*cos(theta)*sin(phi)*v_theta + l*sin(theta)*cos(phi)*v_phi; 
            %          l*sin(theta)*v_theta];

            error('Not implemented.')
            
        end

        function Dq_T = kinetic_energy_gradient_from_velocity(self, q, v)
            
            error('Not implemented.')

        end

        function Dq_T = kinetic_energy_gradient_from_velocity_cyclic(self, x, v)
            
            % theta = x(1);
            % v_phi = v(2);
            % m =  self.MASS;
            % l = self.GEOM(1);
            % 
            % Dq_T = [m*l^2*sin(theta)*cos(theta)*v_phi^2];
            error('Not implemented.')

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(self, q, p)
            
            theta = q(1);
            p_psi = p(2);
            p_phi = p(3);
            J_1 = self.GEOM(1);
            J_3 = self.GEOM(3);

            D_theta_T = p_psi*(p_phi/(J_1*sin(theta))+2*p_phi*cos(theta)^2/(J_1*sin(theta)^3) + p_psi*(J_1-J_3)*cos(theta)/(J_1*J_3*sin(theta)) - ...
                p_psi*cos(theta)*(J_3*cos(theta)^2 + J_1*sin(theta)^2)/(J_1*J_3*sin(theta)^3)) ...
                -p_phi^2*cos(theta)/(J_1*sin(theta)^3);

            Dq_T = [D_theta_T; 0;0];

        end

        function Dq_T = kinetic_energy_gradient_from_momentum_cyclic(self, x, p)
            
            theta = x(1);
            p_psi = p(2);
            p_phi = p(3);
            J_1 = self.GEOM(1);
            J_3 = self.GEOM(3);

            D_theta_T = p_psi*(p_phi/(J_1*sin(theta))+2*p_phi*cos(theta)^2/(J_1*sin(theta)^3) + p_psi*(J_1-J_3)*cos(theta)/(J_1*J_3*sin(theta)) - ...
                p_psi*cos(theta)*(J_3*cos(theta)^2 + J_1*sin(theta)^2)/(J_1*J_3*sin(theta)^3)) ...
                -p_phi^2*cos(theta)/(J_1*sin(theta)^3);

            Dq_T = D_theta_T;

        end

        function L = get_cartesian_angular_momentum_from_momentum(~, q, p)
            
            theta = q(1);
            phi = q(2);
            p_theta = p(1);
            p_phi = p(2);

            L = [sin(phi)*p_theta + cos(theta)/sin(theta)*cos(phi)*p_phi;
                 -p_theta*cos(phi)+ cos(theta)/sin(theta)*sin(phi)*p_phi;
                 p_phi];

        end

        %% Potential functions

        function V_ext = external_potential(self, q)
            % External potential
            theta = q(1);
            m =  self.MASS;
            g =  self.EXT_ACC;
            l = self.GEOM(4);
            V_ext = m*g*l*cos(theta);

        end

        function V_ext = external_potential_cyclic(self, x)
            % External potential
            theta = x(1);
            m =  self.MASS;
            g =  self.EXT_ACC;
            l = self.GEOM(4);
            V_ext = m*g*l*cos(theta);

        end

        function DV_ext = external_potential_gradient(self, q)
            theta = q(1);
            m =  self.MASS;
            g =  self.EXT_ACC;
            l = self.GEOM(1);
            DV_ext = [-m*g*l*sin(theta); 0; 0];
        end

        function DV_ext = external_potential_gradient_cyclic(self, x)
            theta = x(1);
            m =  self.MASS;
            g =  self.EXT_ACC;
            l = self.GEOM(4);
            DV_ext = -m*g*l*sin(theta);
        end

        function D2V_ext = external_potential_hessian(self, q)
            theta = q(1);
            m =  self.MASS;
            g =  self.EXT_ACC;
            l = self.GEOM(4);
            D2V_ext = [-m*g*l*cos(theta), 0, 0; 0, 0, 0; 0,0,0];

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


        %% Constraint on position level

        function g = constraint(~, ~)
           
            g=NaN;

        end

        function Dg = constraint_gradient(~, ~)
           
            Dg=NaN;

        end

        function D2g = constraint_hessian(~, ~, ~)

             D2g=NaN;

        end

        %% Invariant formulations
        %  e.g. for EMS

        % Invariant of the internal potential
        function pi = potential_invariant(~, q, i)

            
                error('system has no invariants for the potential.');
     

        end

            % gradient of potential invariants w.r.t. q
        function DpiDq = potential_invariant_gradient(~, q, i)
            
                error('system has no invariants for the potential.');

        end

        % gradient of potential invariants w.r.t. q
        function D2piDq = potential_invariant_hessian(~, ~, i)

            
                error('system has no invariants for the potential.');
            
        end

        % internal potential computed with the invariant
        function Vs = potential_from_invariant(self, pi, i)
            error('system has no invariants for the potential.');
        end

        % gradient of internal potential w.r.t. the invariant
        function DVsDpi = potential_gradient_from_invariant(self, pi, i)
            error('system has no invariants for the potential.');
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
    function zeta = constraint_invariant(~, ~, ~)

       error('not available.')
    end

    % gradient of the invariant of the position constraint w.r.t. q
    function DzetaDq = constraint_invariant_gradient(~, ~, ~)

        error('not available.')
    end

    % gradient of the invariant of the position constraint w.r.t. q
    function D2zetaDq2 = constraint_invariant_hessian(~, ~, ~)

      error('not available.')
    end

    % position constrained computed with its invariant
      function gs = constraint_from_invariant(~, ~, ~)

            error('not available.')
      end

    % gradient of position constrained w.r.t. its invariant
    function gs = constraint_gradient_from_invariant(~, ~, ~)

        error('not available.')
    end

     function analyzed_quantity = hconvergence_set(~, this_simulation)

        if strcmp(this_simulation.CONV_QUANTITY,'q')
            analyzed_quantity = this_simulation.z(end, 2); %position of 4th particle
        elseif strcmp(this_simulation.CONV_QUANTITY,'p')
            analyzed_quantity = this_simulation.z(end, 4); %momentum of 4th particle
        elseif strcmp(this_simulation.CONV_QUANTITY,'lambda')
            error('not available.')
        else
            error('quantity not yet implemented for convergence analysis.')
        end

    end


    function reference_solution = hconvergence_reference(~, ~, analyzed_quantity)

        reference_solution = analyzed_quantity(:, end, end); %position
      
    end

        %% Animation method
            function give_animation(self, fig, this_simulation)

            DIM = self.DIM;
            q = this_simulation.z(:, 1:2);
            l = self.GEOM(1);
            NT = size(q, 1);

            axis equal
            axis([-1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l]);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            grid on;
            
            ra = self.get_cartesian_coordinates(q(1,:));
            xa = ra(1);
            ya = ra(2);
            za = ra(3);

            for j = 1:NT

                cla(fig);
                hold on

                %% Current position
                x_all = xa;
                y_all = ya;
                z_all = za;
                for k=1:j
                    r = self.get_cartesian_coordinates(q(k,:));
                    x_all = [x_all; r(1)];
                    y_all = [y_all; r(2)];
                    z_all = [z_all; r(3)];
                end
                r_now = self.get_cartesian_coordinates(q(j,:));
                x = r_now(1);
                y = r_now(2);
                z = r_now(3);

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
                plot3(x_all(1:j), y_all(1:j), z_all(1:j), 'k');
                plot3(x, y, z, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
                grid on

                %% current position of the constraint
                x3 = [0; x];
                y3 = [0; y];
                z3 = [0; z];
                view(136, 23)
                drawnow

                end

        end

    end

end
