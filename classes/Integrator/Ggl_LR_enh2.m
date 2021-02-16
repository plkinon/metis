classdef Ggl_LR_enh2 < Integrator
%% Variational integration scheme for GGL-like constrained DAE
%
% - based on constraint on position and velocity level 
%   (GGL-stabilisation)
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - taken from Leimkuhler & Reich , p.189
%
% - constraints are enforced at t_{n+1}
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        
        function self = Ggl_LR_enh2(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 2*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.hasPARA = true;
            self.PARA  = this_simulation.INT_PARA(1);
            self.NAME  = 'GGL-LR-enh2';
        end
        
        function z0 = set_initial_condition(self,this_simulation,this_system)
            
            n  = this_simulation.DIM;
            z0 = [this_simulation.Q_0', (this_system.MASS_MAT * this_simulation.V_0)', self.LM0'];
            
        end
        
        function [resi,tang] = compute_resi_tang(self,zn1,zn,this_problem)
            
            %% Abbreviations
            M  = this_problem.MASS_MAT;
            IM = M\eye(size(M));
            h  = self.DT;
            n  = this_problem.nDOF;
            m  = this_problem.mCONSTRAINTS;
            
            %% Unknows which will be iterated
            qn1     = zn1(1:n);
            pn1     = zn1(n+1:2*n);
            gamman = zn(2*n+1:2*n+m);
            gamma_n1 = zn1(2*n+1:2*n+m);
            kappa   = zn1(2*n+m+1:end);
            G_n1    = this_problem.constraint_gradient(qn1);
            g_n1    = this_problem.constraint(qn1);
            
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            g_n    = this_problem.constraint(qn);
            G_n    = this_problem.constraint_gradient(qn);
            
            %% Alpha-eval. quantities
            al     = self.PARA;
            q_nal  = al*qn1+(1-al)*qn;
            p_bar_n = pn -h/2*G_n'*gamman;
            p_bar_n1 = pn1 - h/2*G_n1'*gamma_n1;
            p_bar_n1mal = (1-al)*p_bar_n1+al*p_bar_n;
            
            G_nal  = this_problem.constraint_gradient(q_nal);
            DV_nal   = this_problem.internal_potential_gradient(q_nal) + this_problem.external_potential_gradient(q_nal);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n    = zeros(n);
            for j = 1:m
                t_n   = t_n + this_problem.constraint_hessian(q_nal,j)*kappa(j);
            end
            
            %% Residual vector 
            resi = [(qn1 - qn - h*IM*p_bar_n1mal - h*IM*G_nal'*kappa);
                    (p_bar_n1 - p_bar_n + h*DV_nal +h*t_n*IM*p_bar_n1mal);
                    (al*g_n1+(1-al)*g_n)                                                      ;
                    G_nal*IM*p_bar_n1mal];

            %% Tangent matrix
            tang=[];
            
        end
        
    end

end