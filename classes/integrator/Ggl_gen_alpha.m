classdef Ggl_gen_alpha < Integrator
%% Midpont-Integration scheme for standard constrained DAE
%
% - based only on constraint on position and velocity level 
%   (GGL-stabilisation)
%
% - for at most quadratic constraint and potential and alpha=1/2 this 
%   method is equivalent to energy-momentum scheme
%
% - for choices alpha={0;1} there are problems
%
% Author: Philipp Kinon
% Date  : 01.12.2020

    properties
        
        ALPHA
        
    end

    methods
        function self = initialise(self,CONFIG,this_problem)
            self.NAME  = 'GGL-Gen-alpha';
            self.nVARS = 2*this_problem.nDOF++2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
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
            gamma_n1  = zn1(2*n+m+1:2*n+2*m);
            G_n1      = this_problem.constraint_gradient(qn1);
            g_n1      = this_problem.constraint(qn1);
            D2g_n1    = this_problem.constraint_hessian(qn1);
            
            %% Known quantities from last time-step
            qn      = zn(1:n);
            pn      = zn(n+1:2*n);
            lambdan = zn(2*n+1:2*n+m);
            gamman  = zn(2*n+m+1:2*n+2*m);
            
            %% Alpha evaluated quantities
            q_nal      = (1-al)*qn + al*qn1;
            q_n1mal    = al*qn + (1-al)*qn1;
            p_n1mal      = al*pn + (1-al)*pn1;
            lambda_nal = (1-al)*lambdan + al*lambda_n1;
            gamma_n1mal = al*gamman + (1-al)*gamma_n1;
            
            DV_nal     = this_problem.potential_gradient(q_nal);
            D2V_nal    = this_problem.potential_hessian(q_nal);
            G_nal      = this_problem.constraint_gradient(q_nal);
            G_n1mal    = this_problem.constraint_gradient(q_n1mal);
            D2g_nal    = this_problem.constraint_hessian(q_nal);
            D2g_n1mal  = this_problem.constraint_hessian(q_n1mal);
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*p_n1mal + h * G_n1mal' *gamma_n1mal ;
                    pn1 - pn + h*DV_nal   + h*lambda_nal*G_nal'       ;
                    g_n1                                              ;
                    G_n1*IM*pn1                                       ];

            %% Tangent matrix
            tang = [eye(n) + h*(1-al)*D2g_n1mal*gamma_n1mal   -h*(1-al)*IM      zeros(n,1)   (1-al)*h*G_n1mal';
                    h*al*D2V_nal + h*al*D2g_nal*lambda_nal    eye(n)            al*h*G_nal'  zeros(n,1)       ; 
                    G_n1                                      zeros(n,1)'       0            0                ;
                    (D2g_n1*IM*pn1)'                          G_n1*IM           0            0               ];
            
        end
        
    end

end