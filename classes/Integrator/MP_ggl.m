classdef MP_ggl < Integrator
%% Gen-alpha-typed-Integration scheme for standard constrained DAE
%
% - based on constraint on position and velocity level 
%   (GGL-stabilisation)
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - not derived from variational principle but simply evaluates RHS of ODE 
%   at t_{n+1/2}  
%
% - constraints are enforced at t_{n+1}
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        
        function self = MP_ggl(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 2*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.NAME  = 'MP-GGL';
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
            gamma_n1  = zn1(2*n+m+1:2*n+2*m);
            G_n1      = this_problem.constraint_gradient(qn1);
            g_n1      = this_problem.constraint(qn1);
          
            %% Known quantities from last time-step
            qn      = zn(1:n);
            pn      = zn(n+1:2*n);
            lambdan = zn(2*n+1:2*n+m);
            gamman  = zn(2*n+m+1:2*n+2*m);
            
            %% Alpha evaluated quantities
            q_n05      = (0.5)*qn + 0.5*qn1;
            p_n05      = 0.5*pn + (0.5)*pn1;
            lambda_n05 = (0.5)*lambdan + 0.5*lambda_n1;
            gamma_n05  = 0.5*gamman + (0.5)*gamma_n1;
            DV_n05     = this_problem.internal_potential_gradient(q_n05) + this_problem.external_potential_gradient(q_n05);
            D2V_n05    = this_problem.internal_potential_hessian(q_n05) + this_problem.external_potential_hessian(q_n05);
            G_n05      = this_problem.constraint_gradient(q_n05);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n1mal    = zeros(n);
            for i = 1:m
                t_n1mal   = t_n1mal + this_problem.constraint_hessian(q_n1mal,i)*gamma_n05(i);
            end
            t_nal = zeros(n);
            for i = 1:m
                t_nal   = t_nal + this_problem.constraint_hessian(q_n05,i)*lambda_n05(i);
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
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*p_n05 + h * G_n05' *gamma_n05     ;
                    pn1 - pn + h*DV_n05   + h*G_nal'*lambda_n05       ;
                    g_n1                                              ;
                    G_n1*IM*pn1                                       ];

            %% Tangent matrix
            tang = [eye(n) + h*(0.5)*t_n1mal                  -h*(0.5)*IM       zeros(n,m)    (0.5)*h*G_n05';
                    h*0.5*D2V_n05 + h*0.5*t_nal               eye(n)            0.5*h*G_nal'  zeros(n,m)       ; 
                    G_n1                                      zeros(n,m)'       zeros(m)      zeros(m)         ;
                    T_n1                                      G_n1*IM           zeros(m)      zeros(m)      ];
            
        end
        
    end

end