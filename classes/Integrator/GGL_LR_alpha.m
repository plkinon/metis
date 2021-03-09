classdef GGL_LR_alpha < Integrator
%% Variational integration scheme for GGL-like constrained DAE
%
% - based on constraint on position and velocity level 
%   (GGL-stabilisation)
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - taken from Leimkuhler & Reich , p.189, using Gen-alpha-method for Psi_h
%
% - not a VI
%
% - constraints are enforced at t_{n+1}
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        
        function self = GGL_LR_alpha(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 2*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.hasPARA = true;
            self.PARA  = this_simulation.INT_PARA(1);
            self.NAME  = 'GGL-LR-alpha';
        end
        
        function z0 = set_initial_condition(self,this_simulation,this_system)
            
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
            lambdan = zn1(2*n+1:2*n+m);
            gamma_n1 = zn1(2*n+m+1:end);
            G_n1    = this_problem.constraint_gradient(qn1);
            g_n1    = this_problem.constraint(qn1);
            
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            
            %% Alpha-eval. quantities
            al     = self.PARA;
            q_nal  = al*qn1+(1-al)*qn;
            G_n    = this_problem.constraint_gradient(qn);
            DV_nal   = this_problem.internal_potential_gradient(q_nal) + this_problem.external_potential_gradient(q_nal);

            %% Residual vector 
            resi = [(qn1 - qn -(1-al)*h*IM*pn1 +(1-al)*(h^2)/2*IM*G_n1'*gamma_n1-al*h*IM*pn+al*(h^2)/2*IM*G_n'*lambdan);
                    (pn1 - pn + h*DV_nal - h/2*G_n1'*gamma_n1 + h/2*G_n'*lambdan);
                    g_n1                                                      ;
                    G_n1*IM*pn1                                              ];

            %% Tangent matrix
            tang=[];
            
        end
        
    end

end