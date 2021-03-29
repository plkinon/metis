classdef RigidBodyMoving < System

    %% Rigid Body moving through space
    methods

        function self = RigidBodyMoving(CONFIG)
            self.nBODIES      = 4;
            self.mCONSTRAINTS = 6;
            self.DIM          = CONFIG.DIM;
            self.MASS         = CONFIG.MASS;
            self.MASS_MAT     = self.MASS*eye(12);
            self.nDOF         = 12;
            self.EXT_ACC      = repmat(CONFIG.EXT_ACC,4,1);
            %self.GEOM(1)      = 0;
            self.nPotentialInvariants  = 0;
            self.nConstraintInvariants = 6;
            self.nVconstraintInvariants = 6;      
                          
        end

        function self = initialise(self, CONFIG, this_integrator)
            % Set initial values
            self.z       = zeros(this_integrator.NT, this_integrator.nVARS);
            self.z(1, :) = [CONFIG.Q_0', (self.MASS_MAT * CONFIG.V_0)', this_integrator.LM0'];
        end
        
        function V_ext = external_potential(self, q)
            V_ext = 0*(self.MASS_MAT*self.EXT_ACC)'*q;
            %% TO DO 
        end
        
        function DV_ext = external_potential_gradient(self,~)
            DV_ext = 0*self.MASS_MAT*self.EXT_ACC;
            %% TO DO
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
            d1 = q(self.DIM+1:2*self.DIM);
            d2 = q(2*self.DIM+1:3*self.DIM);
            d3 = q(3*self.DIM+1:4*self.DIM);
            g1 = 0.5 * (d1' * d1 - 1);
            g2 = 0.5 * (d2' * d2 - 1);
            g3 = 0.5 * (d3' * d3 - 1);
            g4 = 0.5 * (d1' * d2);
            g5 = 0.5 * (d1' * d3);
            g6 = 0.5 * (d2' * d3);
            g  = [g1 ; g2; g3; g4; g5; g6];
        end

        function Dg = constraint_gradient(self, q)
            % Gradient of constraint w.r.t q
            d1 = q(self.DIM+1:2*self.DIM);
            d2 = q(2*self.DIM+1:3*self.DIM);
            d3 = q(3*self.DIM+1:4*self.DIM);
            null = zeros(1,3);
            Dg = [null  d1'   null  null;
                  null  null  d2'   null;
                  null  null  null  d3';
                  null  d2'   d1'   null;
                  null  d3'   null  d1';
                  null  null  d3'   d2'];
        end

        function D2g = constraint_hessian(~,~,m)
            % Hessian of g_1 w.r.t. q
            D2g = zeros(12,12);
            if m == 1
                D2g(4,4) = 1;
                D2g(5,5) = 1;
                D2g(6,6) = 1;
            elseif m==2
                D2g(7,7) =1;
                D2g(8,8) = 1;
                D2g(9,9) = 1;
            elseif m==3
                D2g(10,10) =1;
                D2g(11,11) = 1;
                D2g(12,12) = 1;
            end
        end
        
        function [] = potential_invariant(~,~,~)
            
            error('Problem has only no invariants for the potential.');

        end
        
        function [] = potential_invariant_gradient(~,~,~)
            
            error('Problem has only no invariants for the potential.');
            
        end
        
        %% Invariant formulations
        %  e.g. for EMS
        
        % invariant of the velocity invariant
        function pi2 = vConstraint_invariant(self,q,p,~)
            
            m = self.MASS;
            pi2 = q'*p/m;
            
        end
        
        % gradient of the invariant of the velocity constraint w.r.t. q
        function Dpi2Dq = vConstraint_invariant_gradient_q(self,~,p,~)
                      
            m = self.MASS;
            Dpi2Dq = p/m;
            
        end
        
        % gradient of the invariant of the velocity constraint w.r.t. p
        function Dpi2Dp = vConstraint_invariant_gradient_p(self,q,~,~)
            
            m = self.MASS;
            Dpi2Dp = q/m;
            
        end
        
        function D2piDqDp = vConstraint_invariant_hessian_qp(self,~,~,~)
            
            m = self.MASS;
            D2piDqDp = 1/m*eye(self.DIM);
            
        end
        
        % velocity constraint computed with its invariant
        function gv = Vconstraint_from_invariant(~,pi2,~)
            
            gv = pi2;
            
        end
        
        function DgvDpi = Vconstraint_gradient_from_invariant(~,~,~)
            DgvDpi = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function zeta = constraint_invariant(~,q,i)
                     
            if i == 1
                zeta = q'*q;
            else
                error('Problem has only 1 invariant for the constraint.');
            end
        end
        
        function DzetaDq = constraint_invariant_gradient(~,q,i)
                        
            if i == 1
                DzetaDq = 2*q';
            else
                error('Problem has only 1 invariant for the constraint.');
            end
        end
        
        % gradient of the invariant of the position constraint w.r.t. q
        function D2zetaDq2 = constraint_invariant_hessian(self,~,i)
                                  
            if i == 1
                D2zetaDq2 = 2*eye(self.DIM); 
            else
                error('Problem has only 1 invariant for the constraint.');
            end
        end
        
        function gs = constraint_from_invariant(self,zeta,i)
            
            if i == 1
                gs = 0.5 * (zeta - self.GEOM(1)^2);
            else
                error('Problem has only 1 invariant for the constraint.');
            end
        end
        
        function gs = constraint_gradient_from_invariant(~,~,i)
            
            if i == 1
                gs = 0.5 ;
            else
                error('Problem has only 1 invariant for the constraint.');
            end
        end  
        
        function give_animation(self,fig,this_simulation)
            
                DIM = self.DIM;
                phi = this_simulation.z(:, 1:DIM);
                d1 = this_simulation.z(:, DIM+1:2*DIM);
                d2 = this_simulation.z(:, 2*DIM+1:3*DIM);
                d3 = this_simulation.z(:, 3*DIM+1:4*DIM);
                NT = size(phi,1);

                axis equal
                minx = min(phi(:,1))-1;
                maxx = max(phi(:,1))+1;
                miny = min(phi(:,2))-1;
                maxy = max(phi(:,2))+1;
                minz = min(phi(:,3))-1;
                maxz = max(phi(:,3))+1;
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
                    x1 = x0+d1(j,1);
                    y1 = y0+d1(j,2);
                    z1 = z0+d1(j,3);
                    x2 = x0+d2(j,1);
                    y2 = y0+d2(j,2);
                    z2 = z0+d2(j,3);
                    x3 = x0+d3(j,1);
                    y3 = y0+d3(j,2);
                    z3 = z0+d3(j,3);


                    %% current position of the mass
                    hold on
                    plot3(x0, y0, z0, 'mo', 'MarkerSize', 20, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    plot3(x1, y1, z1, 'mo', 'MarkerSize', 4, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    plot3(x2, y2, z2, 'mo', 'MarkerSize', 4, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    plot3(x3, y3, z3, 'mo', 'MarkerSize', 4, 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [0.75, 0, 0]);
                    
                    grid on

                    %% current position of the constraint
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
                    plot3(xx, yy, zz, 'k', 'linewidth', 1);
                    plot3(xxx, yyy, zzz, 'k', 'linewidth', 1);
                    plot3(xxxx, yyyy, zzzz, 'k', 'linewidth', 1);
                    view(136, 23)
                    
                    drawnow
                    
                end
                
        end

    end

end
