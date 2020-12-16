classdef DoublePendulum < System

    %% Double Pendulum system in 2 or 3 dimensions

    methods

        function self = DoublePendulum(CONFIG)

            self.mCONSTRAINTS = 2;
            self.nBODIES      = 2;
            self.DIM          = CONFIG.DIM;
            self.MASS         = CONFIG.MASS;
            self.nDOF         = self.nBODIES*CONFIG.DIM;
            self.MASS_MAT     = diag([repmat(self.MASS(1),self.DIM,1);repmat(self.MASS(2),self.DIM,1)]);
            self.EXT_ACC      = repmat(CONFIG.EXT_ACC,self.nBODIES,1);
            self.GEOM(1)      = norm(CONFIG.Q_0(1:CONFIG.DIM));   %length of 1st rod
            self.GEOM(2)      = norm(CONFIG.Q_0(CONFIG.DIM+1:2*CONFIG.DIM)-CONFIG.Q_0(1:CONFIG.DIM)); %length of 2nd rod
            self.nPotentialInvariants  = 0;
            self.nConstraintInvariants = 2;
        end

        function self = initialise(self, CONFIG, this_integrator)
            % Set initial values
            self.z       = zeros(this_integrator.NT, this_integrator.nVARS);
            self.z(1, :) = [CONFIG.Q_0', (self.MASS_MAT * CONFIG.V_0)', this_integrator.LM0'];
        end
        
        function V_ext = external_potential(self, q)
            V_ext = (self.MASS_MAT*self.EXT_ACC)'*q;
        end
        
        function DV_ext = external_potential_gradient(self,~)
            DV_ext = self.MASS_MAT*self.EXT_ACC;
        end
        
        function D2V_ext = external_potential_hessian(~,q)
            D2V_ext = zeros(size(q,1));
        end

        function V_int = internal_potential( ~, ~)
            V_int = 0;
        end
        
        function DV_int = internal_potential_gradient(~,q)
            DV_int = zeros(size(q));
        end
        
        function D2V_int = internal_potential_hessian(~,q)
            D2V_int = zeros(size(q,1));
        end

        function g = constraint(self, q)
            % Constraint on position level
            q1 = q(1:self.DIM);
            q2 = q(self.DIM+1:2*self.DIM);
            g1 = 0.5 * (q1' * q1 - self.GEOM(1)^2);
            g2 = 0.5 * ((q2-q1)' * (q2-q1) -self.GEOM(2)^2);
            g  = [g1 ; g2];
            
        end

        function Dg = constraint_gradient(self, q)
            % Gradient of constraint w.r.t q
            q1 = q(1:self.DIM);
            q2 = q(self.DIM+1:2*self.DIM);
            Dg = [q1' zeros(size(q1))'; -(q2-q1)' (q2-q1)'];
            
        end
        
        function D2g = constraint_hessian(self,~,m)
            
            if m == 1
                % Hessian of g_1 w.r.t. q
                D2g = [eye(self.DIM)  zeros(self.DIM); zeros(self.DIM) zeros(self.DIM)]; 
            elseif 2
                % Hessian of g_2 w.r.t. q
                D2g = [eye(self.DIM) -eye(self.DIM); -eye(self.DIM) eye(self.DIM)];
            end
            
        end
        
        function [] = potential_invariant(~,~,~)
            
            error('Problem has only no invariants for the potential.');

        end
        
        function [] = potential_invariant_gradient(~,~,~)
            
            error('Problem has only no invariants for the potential.');
            
        end
        
        function zeta = constraint_invariant(self,q,i)
            
            q1 = q(1:self.DIM);
            q2 = q(self.DIM+1:2*self.DIM);
            
            if i == 1
                zeta = q1'*q1;
            elseif i == 2
                zeta = (q2-q1)'*(q2-q1);
            else
                error('Problem has only 2 invariants for the constraint.');
            end
        end
        
        function DzetaDq = constraint_invariant_gradient(self,q,i)
            
            q1 = q(1:self.DIM);
            q2 = q(self.DIM+1:2*self.DIM);
            
            if i == 1
                DzetaDq = [2*q1; zeros(self.DIM,1)];
            elseif i == 2
                DzetaDq = [-2*(q2-q1); 2*(q2-q1)];
            else
                error('Problem has only 2 invariants for the constraint.');
            end
        end
        
        function give_animation(self,fig)
            
            
                DIM = self.DIM;
                q1 = self.z(:, 1:DIM);
                q2 = self.z(:,DIM+1:2*DIM);
                l = self.GEOM(1)+self.GEOM(2);
                NT = size(q1,1);

                axis equal
                axis([-1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l, -1.1 * l, 1.1 * l]);
                xlabel('x');
                ylabel('y');
                zlabel('z');
                grid on;

                xa1 = q1(1, 1);
                xa2 = q2(1, 1);
                ya1 = q1(1, 2);
                ya2 = q2(1, 2);
                if DIM == 3
                    za1 = q1(1, 3);
                    za2 = q2(1, 3);
                else
                    za1 = 0;
                    za2 = 0;
                end

                for j = 1:NT

                    cla(fig);
                    hold on

                    %% Current position
                    x1 = q1(j, 1);
                    x2 = q2(j, 1);
                    y1 = q1(j, 2);
                    y2 = q2(j, 2);
                    
                    if DIM == 3
                        z1 = q1(j, 3);
                        z2 = q2(j, 3);
                    else
                        z1 = 0;
                        z2 = 0;
                    end

                    %% Reference sphere
                    plot3(xa1, ya1, za1, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    hold on
                    plot3(xa2, ya2, za2, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [0.75, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    hold on

                    %% Reference constraint
                    x3 = [0; xa1; xa2];
                    y3 = [0; ya1; ya2];
                    z3 = [0; za1; za2];
                    plot3(x3, y3, z3, 'k', 'LineWidth', 1);

                    %% current position of the mass
                    hold on
                    if DIM == 3
                        plot3(q1(1:j, 1), q1(1:j, 2), q1(1:j, 3), 'k');
                        plot3(q2(1:j, 1), q2(1:j, 2), q2(1:j, 3), 'k');
                    else
                        plot3(q1(1:j, 1), q1(1:j, 2), zeros(j, 1), 'k');
                        plot3(q2(1:j, 1), q2(1:j, 2), zeros(j, 1), 'k');
                    end
                    plot3(x1, y1, z1, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    plot3(x2, y2, z2, 'mo', 'MarkerSize', 10, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    grid on

                    %% current position of the constraint
                    x3 = [0; x1; x2];
                    y3 = [0; y1; y2];
                    z3 = [0; z1; z2];
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
