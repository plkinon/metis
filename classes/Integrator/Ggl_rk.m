classdef Ggl_rk < Integrator
%% Runge-Kutta typed scheme for GGL-like constrained DAE
%
% - based on constraint on position and velocity level 
%   (GGL-stabilisation)
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - derived from variational principle by Peter Betsch
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        
        function self = Ggl_rk(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 3*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.PARA  = [1 1]; %[The:round theta  theta: vartheta];
            self.NAME  = 'GGL-RK';
        end
        
        function z0 = set_initial_condition(self,this_simulation,this_system)
            
            M  = this_system.MASS_MAT;
            v0 = this_simulation.V_0;
            q0 = this_simulation.Q_0;
            n  = this_system.nDOF;
            m  = this_system.mCONSTRAINTS;
            h  = self.DT;
            The     = self.PARA(1);
            lambda0 = self.LM0(1:m);
            gamma0 = self.LM0(m+1:2*m);
            G_0    = this_system.constraint_gradient(q0);
            DV_0   = this_system.internal_potential_gradient(q0) + this_system.external_potential_gradient(q0);
            t_0    = zeros(n);
            for j = 1:m
                t_0   = t_0 + this_system.constraint_hessian(q0,j)*gamma0(j);
            end
            p0 = (M+h*(1-The)*t_0)*v0 + h*((1-The)*DV_0 + (1-theta)*G_0'*lambda0); 
            %p0 = M*v0;
            z0 = [q0', p0', v0' , self.LM0'];
            
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
            vn1     = zn1(2*n+1:3*n);
            lambdan = zn1(3*n+1:3*n+m);
            gamman  = zn1(3*n+m+1:end);
            G_n1    = this_problem.constraint_gradient(qn1);
            g_n1    = this_problem.constraint(qn1);
            
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
            G_n    = this_problem.constraint_gradient(qn);
            g_n    = this_problem.constraint(qn);
            
            %% Quantities at t_n+theta
            theta   = self.PARA(2);
            The     = self.PARA(1);
            q_nt    = (1-The)*qn + The*qn1;
            DV_nt   = this_problem.internal_potential_gradient(q_nt) + this_problem.external_potential_gradient(q_nt);
            G_nt    = this_problem.constraint_gradient(q_nt);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n    = zeros(n);
            for j = 1:m
                t_n   = t_n + this_problem.constraint_hessian(q_nt,j)*gamman(j);
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*vn1 - h*IM*G_nt'*gamman                  ;
                    pn1 - pn + h*DV_nt + h*((1-theta)*G_n'+theta*G_n1')*lambdan + h*t_n*vn1;
                    (M + h * (1-The) * t_n) * vn1 - pn + h * ( (1-The) * DV_nt + (1-theta) * G_n' * lambdan) ;
                    theta*g_n1 + (1-theta)*g_n                                             ;
                    G_nt*vn1                                              ];

            %% Tangent matrix
            tang = [];
%             tang = [eye(n) - h*IM*t_n1          -h*IM         zeros(n,m)   -h*IM*G_n1' ;
%                     zeros(n)                    eye(n)       h*G_n'        zeros(n,m)  ; 
%                     G_n1                        zeros(n,m)'  zeros(m)      zeros(m)    ;
%                     T_n1                        G_n1*IM      zeros(m)      zeros(m)    ];
%             
        end
        
    end

end