classdef Ggl_std < Integrator
%% GGL integration scheme for constrained DAE
%
% - based constraints on position and velocity level
%
% - standard stabilisation scheme by Gear, Gupta & Leimkuhler where
%   lambda_n serves as an unknown and therefore it differs from GGL with
%   Gen-alpha=0
%
% Author: Philipp Kinon
% Date  : 01.12.2020
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
            
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            G_n    = this_problem.constraint_gradient(qn);
            D2g    = this_problem.constraint_hessian(qn);
            DV_n   = this_problem.potential_gradient(qn);
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*pn1 + h*gamman1*G_n1' ;
                    pn1 - pn + h*DV_n   + h*lambdan*G_n'  ;
                    g_n1                                  ;
                    G_n1*IM*pn1                          ];

            %% Tangent matrix
            tang = [eye(n) + h*D2g*gamman1      -h*IM          zeros(n,1)    h*G_n1'     ;
                    zeros(n)                    eye(n)         h*G_n'        zeros(n,1)  ; 
                    G_n1                        zeros(n,1)'    0             0           ;
                    (D2g*IM*pn1)'               G_n1*IM        0             0           ];
            
        end
        
    end

end