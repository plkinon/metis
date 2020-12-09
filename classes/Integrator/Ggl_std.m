classdef Ggl_std < Integrator
%% GGL integration scheme for constrained DAE
%
% - based constraints on position and velocity level
%
% - independent momenta variables (Hamilton-Potryagin approach)
%
% - standard stabilisation scheme by Gear, Gupta & Leimkuhler where
%   lambda_n serves as an unknown
%
% - not derived from variational principle: constraints evaluated at
%   t_{n+1} and ODE-RHS for q at t_{n+1} and for p at t_n
%
% - bad performance!
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        function self = initialise(self,~,this_problem)
            self.NAME  = 'GGL-std ';
            self.nVARS = 2*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
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
            gamman1 = zn1(2*n+m+1:end);
            G_n1    = this_problem.constraint_gradient(qn1);
            g_n1    = this_problem.constraint(qn1);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n1    = zeros(n);
            for i = 1:m
                t_n1   = t_n1 + this_problem.constraint_hessian(qn1,i)*gamman1(i);
            end
            
            % Hessian of constraints are multiplied by inverse MassMat and
            % pn1 for each constraint to avoid 3rd order tensors
            T_n1 = zeros(m,n);
            for i = 1:m
                tmp = this_problem.constraint_hessian(qn1,i);
                for k = 1:n
                    T_n1(i,k) = tmp(:,k)'*IM*pn1;
                end
            end
            
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            G_n    = this_problem.constraint_gradient(qn);
            DV_n   = this_problem.potential_gradient(qn);
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*pn1 + h*G_n1'*gamman1 ;
                    pn1 - pn + h*DV_n   + h*G_n'*lambdan  ;
                    g_n1                                  ;
                    G_n1*IM*pn1                          ];

            %% Tangent matrix
            tang = [eye(n) + h*t_n1             -h*IM          zeros(n,m)    h*G_n1'     ;
                    zeros(n)                    eye(n)         h*G_n'        zeros(n,m)  ; 
                    G_n1                        zeros(n,m)'    zeros(m)      zeros(m)           ;
                    T_n1                        G_n1*IM        zeros(m)      zeros(m)    ];
            
        end
        
    end

end