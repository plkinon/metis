classdef GGL_VI < Integrator
%% Variational integration scheme for GGL-like constrained DAE
%
% - based on constraint on position and momentum level 
%
% - independent momenta variables (Livens approach)
%
% - derived from variational principle by Peter Betsch (1st attempt for
%   new 'GGL-functional'
%
% - not symplectic
%
% - constraints are enforced at t_{n+1}
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        
        function self = GGL_VI(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 2*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.INDI_VELO = false;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.hasPARA = false;
            self.NAME  = 'GGL-VI (1st attempt) ';
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
            for l = 1:m
                tmp = this_problem.constraint_hessian(qn1,l);
                for k = 1:n
                    T_n1(l,k) = tmp(:,k)'*IM*pn1;
                end
            end
            
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            gamman = zn(2*n+m+1:end);
            G_n    = this_problem.constraint_gradient(qn);
            DV_n   = this_problem.internal_potential_gradient(qn) + this_problem.external_potential_gradient(qn);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n    = zeros(n);
            for j = 1:m
                t_n   = t_n + this_problem.constraint_hessian(qn,j)*gamman(j);
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*pn1 - h*IM*G_n1'*gamman1                  ;
                    pn1 - pn + h*DV_n + h*G_n'*lambdan + h*t_n*IM*pn          ;
                    g_n1                                                      ;
                    G_n1*IM*pn1                                              ];

            %% Tangent matrix
            tang = [eye(n) - h*IM*t_n1          -h*IM         zeros(n,m)   -h*IM*G_n1' ;
                    zeros(n)                    eye(n)       h*G_n'        zeros(n,m)  ; 
                    G_n1                        zeros(n,m)'  zeros(m)      zeros(m)    ;
                    T_n1                        G_n1*IM      zeros(m)      zeros(m)    ];
            
        end
        
    end

end