classdef Pendulum < System

    %% Pendulum system in 2 or 3 dimensions

    methods

        function self = Pendulum(CONFIG)

            self.MASS_MAT     = eye(CONFIG.DIM);
            self.DIM          = CONFIG.DIM;
            self.nDOF         = CONFIG.nDOF;
            self.mCONSTRAINTS = 1;
            self.EXT_ACC      = reshape(CONFIG.EXT_ACC,[CONFIG.nDOF,1]);
            self.GEOM(1)      = 1; %length of pendulum
        end

        function self = initialise(self, CONFIG, this_integrator)
            % Set initial values
            self.z       = zeros(this_integrator.NT, this_integrator.nVARS);
            self.z(1, :) = [CONFIG.Q_0', (self.MASS_MAT * CONFIG.V_0)', this_integrator.LM0'];
        end
        
        function V = potential(self, q)
            V = (self.MASS_MAT*self.EXT_ACC)'*q;
        end
        
        function DV = potential_gradient(self,~)
            DV = self.MASS_MAT*self.EXT_ACC;
        end
        
        function D2V = potential_hessian(~,q)
            D2V = zeros(size(q,1));
        end

        function g = constraint(self, q)
            % Constraint on position level
            g = 0.5 * (q' * q - self.GEOM(1)^2);
        end

        function Dg = constraint_gradient(~, q)
            % Gradient of constraint w.r.t q
            Dg = q';
        end

        function D2g = constraint_hessian(~,q)
            % Hessian of g w.r.t. q
            D2g = eye(size(q,1));
        end

    end

end
