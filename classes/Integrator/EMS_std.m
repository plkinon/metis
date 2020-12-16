classdef EMS_std < Integrator
%% Energy_Momentum-Integration scheme for standard constrained DAE
%
% - based only on constraint on position level
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - not derived from variational principle
% 
% - taken from Gonzales 1999
%
% - uses standard gradient for ext. potential and discrete gradient for
%   internal potential and constraint
%
% Author: Philipp Kinon
% Date  : 16.12.2020

    methods
        function self = initialise(self,~,this_problem)
            self.NAME  = 'EMS-std';
            self.nVARS = 2*this_problem.nDOF+this_problem.mCONSTRAINTS;
            self.LM0   = zeros(this_problem.mCONSTRAINTS,1);
        end
        
        function [resi,tang] = compute_resi_tang(self,zn1,zn,this_problem)
            
            %% Abbreviations
            M  = this_problem.MASS_MAT;
            IM = M\eye(size(M));
            h  = self.DT;
            n  = this_problem.nDOF;
            m  = this_problem.mCONSTRAINTS;
            
            %% Unknows which will be iterated
            qn1       = zn1(1:n);
            pn1       = zn1(n+1:2*n);
            lambda_n1 = zn1(2*n+1:2*n+m);
            G_n1      = this_problem.constraint_gradient(qn1);
            g_n1      = this_problem.constraint(qn1);
            
            %% Known quantities from last time-step
            qn      = zn(1:n);
            pn      = zn(n+1:2*n);
            lambdan = zn(2*n+1:2*n+m);
            
            %% MP evaluated quantities
            q_n05      = 0.5*(qn + qn1);
            p_n05      = 0.5*(pn + pn1);
            lambda_n05 = 0.5*(lambdan + lambda_n1);
            DVext_n05  = this_problem.external_potential_gradient(q_n05);
            D2Vext_n05 = this_problem.external_potential_hessian(q_n05);
            G_n05      = this_problem.constraint_gradient(q_n05);
            t_n05      = zeros(n);
            
            %% Discrete gradients
            DG_Vint = zeros(n,1);
            for i = 1:this_problem.nPotentialInvariants
                pi_n     = this_problem.potential_invariant(qn,i);
                pi_n1    = this_problem.potential_invariant(qn1,i);
                Vs_n     = this_problem.potential_from_invariant(pi_n,i);
                Vs_n1    = this_problem.potential_from_invariant(pi_n1,i);
                DPiq_n05 = this_problem.potential_invariant_gradient(q_n05,i);
                DG_Vint  = DG_Vint + (Vs_n1-Vs_n)/(pi_n1-pi_n)*DPiq_n05;
            end
            
            DG_g = zeros(m,n);
            for j = 1:this_problem.nConstraintInvariants
                zeta_n = this_problem.constraint_invariant(qn,j);
                zeta_n1 = this_problem.constraint_invariant(qn1,j);
                gs_n = this_problem.constraint_from_invariant(zeta_n,j);
                gs_n1 = this_problem.constraint_from_invariant(zeta_n1,j);
                DzetaDq_n05 = this_problem.constraint_invariant_gradient(q_n05,j);
                DG_g(j,:) = (gs_n1 - gs_n)/(zeta_n1 - zeta_n)*DzetaDq_n05';
            end
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            for i = 1:m
                t_n05   = t_n05 + this_problem.constraint_hessian(q_n05,m)*lambda_n1(m);
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*p_n05                                     ;
                    pn1 - pn + h*DVext_n05   + h*DG_Vint + h*DG_g'*lambda_n1  ;
                    g_n1                                                     ];

            %% Tangent matrix
            tang = [eye(n)                                    -h*0.5*IM      zeros(n,m) ;
                    h*0.5*D2Vext_n05 + h*0.5*t_n05                 eye(n)            h*G_n05'; 
                    G_n1                                      zeros(n,m)'       zeros(m)          ];
            
        end
        
    end

end