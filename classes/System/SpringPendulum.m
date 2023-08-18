%% Class: SpringPendulum
classdef SpringPendulum < System
    % Spring pendulum with nonlinear elastic spring law.

    properties

    end

    methods

        function self = SpringPendulum(CONFIG)

            self.mCONSTRAINTS = 0;
            self.nBODIES = 1;
            self.DIM = CONFIG.DIM;
            self.MASS = CONFIG.MASS;
            self.nDOF = 3;
            self.mMixedQuantities = 1;
            self.MASS_MAT = eye(3)*self.MASS;
            self.EXT_ACC = CONFIG.EXT_ACC;

           % Resting lengths of the springs
            self.GEOM(1) = 1;
            %k (stiffness)
            self.GEOM(2) = 10000;

            self.nPotentialInvariants = 1;

            self.DISS_MAT = zeros(3,3);
            self.isCyclicCoordinate = [false;false;false];

        end
        
        function M = get_mass_matrix(self, ~)
            
            M = self.MASS_MAT;

        end

        function Dq_T = kinetic_energy_gradient_from_velocity(~, ~, ~)
            
            Dq_T = zeros(3,1);

        end

        function Dq_T = kinetic_energy_gradient_from_momentum(~, ~, ~)
            
            Dq_T = zeros(3,1);

        end

        %% Potential functions

        function V_ext = external_potential(self, q)
            % External potential

            m =  self.MASS;
            b =  self.EXT_ACC;
            V_ext = -m*b'*q;

        end

        function DV_ext = external_potential_gradient(self, ~)

            m =  self.MASS;
            b =  self.EXT_ACC;
            DV_ext = -m*b;

        end

        function D2V_ext = external_potential_hessian(~, ~)
            
            D2V_ext = zeros(3,3);

        end

        function V_int = internal_potential(self, q)
            
            l0 = self.GEOM(1);
            k = self.GEOM(2);

            C = (q'*q)/(l0^2);
            V_int = k*l0/4*(C-log(C)-1);

        end

        function V_int = internal_potential_from_mixed_quantity(self, C)
            
            l0 = self.GEOM(1);
            k = self.GEOM(2);

            V_int = k*l0/4*(C-log(C)-1);
            
        end

        function DV_int = internal_potential_gradient(self, q)
            l0 = self.GEOM(1);
            k = self.GEOM(2);

            eps = (q'*q-l0^2)/(2*l0^2);

            DV_int = k*l0/2*(1- 1/(2*eps+1)) * q/l0^2 ;
        end

        function DV_int = internal_potential_gradient_from_mixed_quantity(self, C)
            l0 = self.GEOM(1);
            k = self.GEOM(2);

            DV_int = k*l0/4*(1- 1/C) ;
        end
        
        function C_q = mixed_quantity(self,q)
            l0 = self.GEOM(1);
            C_q = (q'*q)/(l0^2);
        end

        function D_C_q = mixed_quantity_gradient(self,q)
            l0 = self.GEOM(1);
            D_C_q = 2*q/(l0^2) ;

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
        function pi = potential_invariant(~, q, ~)

            pi = NaN;

        end

            % gradient of potential invariants w.r.t. q
        function DpiDq = potential_invariant_gradient(~, q, ~)

            DpiDq = NaN;
        end

        % gradient of potential invariants w.r.t. q
        function D2piDq = potential_invariant_hessian(~, ~, ~)

                 D2piDq = NaN;
            
        end

        % internal potential computed with the invariant
        function Vs = potential_from_invariant(self, pi, ~)

            Vs = NaN;
        end

        % gradient of internal potential w.r.t. the invariant
        function DVsDpi = potential_gradient_from_invariant(self, pi, ~)
          
            DVsDpi = NaN;
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
                                                                                                q = this_simulation.z(:, 1:DIM);
                                                                                                l = 2* self.GEOM(1);
                                                                                                NT = size(q, 1);

                                                                                                axis equal
                                                                                                xmin = min(min(q(:, 1)));
                                                                                                ymin = min(min(q(:, 2)));
                                                                                                zmin = min(min(q(:, 3)));
                                                                                                xmax = max(max(q(:, 1)));
                                                                                                ymax = max(max(q(:, 2)));
                                                                                                zmax = max(max(q(:, 3)));
                                                                                                axis([1.1 * xmin, 1.1 * xmax, 1.1 * ymin, 1.1 * ymax, 1.1 * zmin, 1.1 * zmax]);
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
                                                                                                    %plot3(x3, y3, z3, 'k', 'LineWidth', 1);
                                                                                                    tmp = gfx_springelement([0,0,0], [xa, ya, za], 1, 0.05, 7);

                                                                                                            %plot3(xx3, yy3, zz3, 'k--');
                                                                                                            plot3(tmp(1, :), tmp(2, :), tmp(3, :), 'k', 'LineWidth', 1.5);
                                                                                                    %% current position of the mass
                                                                                                    hold on
                                                                                                    if DIM == 3
                                                                                                        plot3(q(1:j, 1), q(1:j, 2), q(1:j, 3), 'color', '#4477AA');
                                                                                                    else
                                                                                                        plot3(q(1:j, 1), q(1:j, 2), zeros(j, 1), 'color', '#4477AA');
                                                                                                    end
                                                                                                    plot3(x, y, z, 'mo', 'MarkerSize', 40, 'MarkerEdgeColor', '#4477AA', 'MarkerFaceColor', '#4477AA');
                                                                                                    grid on

                                                                                                    %% current position of the constraint
                                                                                                    x3 = [0; x];
                                                                                                    y3 = [0; y];
                                                                                                    z3 = [0; z];
                                                                                                    %plot3(x3, y3, z3, 'k', 'linewidth', 1);

%% Current position of the springs


                                                                                                            tmp = gfx_springelement([0,0,0], [x, y, z], 1, 0.05, 7);

                                                                                                            %plot3(xx3, yy3, zz3, 'k--');
                                                                                                            plot3(tmp(1, :), tmp(2, :), tmp(3, :), 'k', 'LineWidth', 1.5);

                                                                                                            %plot3(xxxx3, yyyy3, zzzz3, 'k--');

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
