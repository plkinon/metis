classdef GGL_VI_mod < Integrator
%% Variational integration scheme for GGL-like constrained DAE
%
% - based on constraint on position and momentum level 
%
% - independent momenta variables (Livens approach)
%
% - derived from variational principle by Peter Betsch (easy (2nd) attempt for
%   new 'GGL-functional'
%
% - symplectic 
%
% - constraints are enforced at t_{n+1}
%
% Author: Philipp Kinon
% Date  : 30.11.2021

    methods
        
        function self = GGL_VI_mod(this_simulation,this_problem)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 3*this_problem.nDOF+2*this_problem.mCONSTRAINTS;
            self.INDI_VELO = true;
            self.LM0   = zeros(2*this_problem.mCONSTRAINTS,1);
            self.hasPARA = false;
            self.NAME  = 'GGL-VI (modified) ';
        end
        
        function z0 = set_initial_condition(self,this_simulation,this_system)
            
            M  = this_system.MASS_MAT;
            z0 = [this_simulation.Q_0', (M*this_simulation.V_0)',  this_simulation.V_0', self.LM0'];
            
        end
        
        function z_rearranged = rearrange_unknowns(~,this_simulation,this_problem)
            
            % v_n is an unknown of this scheme, has to be shifted backwards
            % by 1 after computation
            n  = this_problem.nDOF;
            z_rearranged = this_simulation.z;
            z_rearranged(1:(end-1),2*n+1:3*n) = this_simulation.z(2:end,2*n+1:3*n);
            z_rearranged(end,2*n+1:3*n) = NaN;
            
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
            vn      = zn1(2*n+1:3*n);
            lambdan = zn1(3*n+1:3*n+m);
            gamman1 = zn1(3*n+m+1:end);
            g_n1    = this_problem.constraint(qn1);
            G_n1  = this_problem.constraint_gradient(qn1);
                      
            %% Known quantities from last time-step
            qn     = zn(1:n);
            pn     = zn(n+1:2*n);
            G_n    = this_problem.constraint_gradient(qn);
            DV_n   = this_problem.internal_potential_gradient(qn) + this_problem.external_potential_gradient(qn);
            q_bar  = qn + h*vn;
            G_bar  = this_problem.constraint_gradient(q_bar);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n_bar    = zeros(n);
            for j = 1:m
                t_n_bar   = t_n_bar + this_problem.constraint_hessian(q_bar,j)*gamman1(j);
            end
            T_bar = zeros(m,n);
            for l = 1:m
                tmp = this_problem.constraint_hessian(q_bar,l);
                for k = 1:n
                    T_bar(l,k) = tmp(:,k)'*IM*pn1;
                end
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*vn - h*IM*G_bar'*gamman1                  ;
                    pn1 - pn + h*DV_n + h*G_n'*lambdan + h*t_n_bar*IM*pn1          ;
                    M*vn - pn1 - h*t_n_bar*IM*pn1;
                    g_n1                                                      ;
                    G_bar*IM*pn1                                              ];

            %% Tangent matrix
            tang = [[eye(n)                  zeros(n)                  -h*eye(n)-h^2*IM*t_n_bar      zeros(n,m)    -h*IM*G_bar'];
                    [zeros(n)                eye(n)+h*t_n_bar*IM       zeros(n)                      h*G_n'        h*T_bar'    ]; 
                    [zeros(n)                -eye(n)-h*t_n_bar*IM      M                             zeros(n,m)    -h*T_bar'   ];
                    [G_n1                    zeros(n,m)'               zeros(n,m)'                   zeros(m)      zeros(m)    ];
                    [zeros(n,m)'             G_bar*IM                  T_bar*h                       zeros(m)      zeros(m)    ]];
            
        end
        
    end

end