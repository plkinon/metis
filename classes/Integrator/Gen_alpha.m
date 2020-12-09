classdef Gen_alpha < Integrator
%% Gen-Alpha-typed-Integration scheme for standard constrained DAE
%
% - based only on constraint on position level
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - not derived from variational principle but simply evaluates RHS at
%   t_{n+alpha} or t_{n+1-alpha} respectively
% 
% - setting alpha apart from 1/2 leads to problems
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    properties
        
        ALPHA
        
    end

    methods
        function self = initialise(self,CONFIG,this_problem)
            self.NAME  = 'Gen-alpha';
            self.nVARS = 2*this_problem.nDOF+this_problem.mCONSTRAINTS;
            self.LM0   = zeros(this_problem.mCONSTRAINTS,1);
            self.ALPHA = CONFIG.INTEGRATION_VARIABLES(1);
        end
        
        function [resi,tang] = compute_resi_tang(self,zn1,zn,this_problem)
            
            %% Abbreviations
            M  = this_problem.MASS_MAT;
            IM = M\eye(size(M));
            h  = self.DT;
            n  = this_problem.nDOF;
            m  = this_problem.mCONSTRAINTS;
            al = self.ALPHA;
            
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
            
            %% Alpha evaluated quantities
            q_nal      = (1-al)*qn + al*qn1;
            p_nal      = al*pn + (1-al)*pn1;
            lambda_nal = (1-al)*lambdan + al*lambda_n1;
            DV_nal     = this_problem.potential_gradient(q_nal);
            D2V_nal    = this_problem.potential_hessian(q_nal);
            G_nal      = this_problem.constraint_gradient(q_nal);
            t_nal    = zeros(n);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            for i = 1:m
                t_nal   = t_nal + this_problem.constraint_hessian(q_nal,m)*lambda_nal(m);
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*p_nal                        ;
                    pn1 - pn + h*DV_nal   + h*G_nal'*lambda_nal  ;
                    g_n1                                         ];

            %% Tangent matrix
            tang = [eye(n)                                    -h*(1-al)*IM      zeros(n,m) ;
                    h*al*D2V_nal + h*al*t_nal                 eye(n)            al*h*G_nal'; 
                    G_n1                                      zeros(n,m)'       zeros(m)          ];
            
        end
        
    end

end