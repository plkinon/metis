classdef Ggl_variational < Integrator
% GGL Variational Integrator without starting procedure (1st attempt by PB)

    methods
        function self = initialise(self,~,this_problem)
            self.NAME  = 'GGL-VI (1st attempt) ';
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
            gamman = zn(2*n+m+1:end);
            G_n    = this_problem.constraint_gradient(qn);
            D2g    = this_problem.constraint_hessian(qn);
            DV_n   = this_problem.potential_gradient(qn);
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*pn1 - h*IM*gamman1*G_n1'                  ;
                    pn1 - pn + h*DV_n + h*lambdan*G_n' + D2g*h*gamman*IM*pn   ;
                    g_n1                                                      ;
                    G_n1*IM*pn1                                              ];

            %% Tangent matrix
            tang = [eye(n) - h*IM*D2g*gamman1   -h*IM        zeros(n,1)    -h*IM*G_n1' ;
                    zeros(n)                    eye(n)       h*G_n'        zeros(n,1)  ; 
                    G_n1                        zeros(n,1)'  0             0           ;
                    (D2g*IM*pn1)'               G_n1*IM      0             0          ];
            
        end
        
    end

end