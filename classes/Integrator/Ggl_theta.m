classdef Ggl_theta < Integrator
%% Runge-Kutta typed scheme for GGL-like constrained DAE
%
% - based on constraint on position and velocity level 
%   (GGL-stabilisation)
%
% - independent momentum variables (Hamilton Potryagin approach)
%
% - derived from variational principle by Peter Betsch
%
% Author: Philipp Kinon
% Date  : 28.01.2021

    methods
        
        function self = Ggl_theta(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 2*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.PARA  = 0.5;
            self.NAME  = 'GGL-theta';
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
            lambda_n1 = zn1(2*n+1:2*n+m);
            gamman  = zn1(2*n+m+1:end);
            G_n1    = this_problem.constraint_gradient(qn1);
            g_n1    = this_problem.constraint(qn1);
            
            % Hessian of constraints are multiplied by inverse MassMat and
            % pn1 for each constraint to avoid 3rd order tensors
            T_n1 = zeros(m,n);
            %t_n1 = zeros(n);
            for l = 1:m
                tmp = this_problem.constraint_hessian(qn1,l);
                %t_n1   = t_n1 + this_problem.constraint_hessian(qn1,l)*lambdan(l);
                for k = 1:n
                    T_n1(l,k) = tmp(:,k)'*IM*pn1;
                end
            end
            
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            lambdan = zn(2*n+1:2*n+m);
            G_n    = this_problem.constraint_gradient(qn);
            
            %% Quantities at t_n+theta
            theta   = self.PARA(1);
            q_nt    = (1-theta)*qn + theta*qn1;
            p_n1mt  = theta*pn + (1-theta)*pn1;
            lambda_nt = (1-theta)*lambdan + theta*lambda_n1;
            g_nt    = this_problem.constraint(q_nt);
            DV_nt   = this_problem.internal_potential_gradient(q_nt) + this_problem.external_potential_gradient(q_nt);
            G_nt    = this_problem.constraint_gradient(q_nt);
            D2V_nt  = this_problem.internal_potential_hessian(q_nt) + this_problem.external_potential_hessian(q_nt);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_nt    = zeros(n);
            T_nt    = zeros(m,n);
            for j = 1:m
                tmp2 = this_problem.constraint_hessian(q_nt,j);
                t_nt   = t_nt + this_problem.constraint_hessian(q_nt,j)*gamman(j);
                for k = 1:n
                    T_nt(j,k) = tmp2(:,k)'*IM*p_n1mt;
                end
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*p_n1mt - h*IM*G_nt'*gamman                  ;
                    pn1 - pn + h*DV_nt + h*G_nt'*lambda_nt + h*t_nt*IM*p_n1mt;
                    g_nt                                              ;
                    G_nt*IM*p_n1mt                                              ];

            %% Tangent matrix
            tang = [];
            %tang = [eye(n) - h*IM*The*t_nt     zeros(n)          -h*eye(n)                   zeros(n,m)                      -h*IM*G_nt' ;
            %        h*D2V_nt*The + theta*t_n1  eye(n)            h*t_nt                      h*((1-theta)*G_n'+theta*G_n1')   h*T_nt'  ;
            %        h*(1-The)*D2V_nt           zeros(n)          (M + h * (1-The) * t_nt)    h*(1-theta) * G_n'                h*(1-The)*T_nt';     
            %        theta*G_n1                 zeros(n,m)'       zeros(n,m)'                 zeros(m)                        zeros(m)    ;
            %s        The*T_nt                   zeros(n,m)'       G_nt          zeros(m)      zeros(m)    ];
             
        end
        
    end

end