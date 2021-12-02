classdef MP_std < Integrator
%% Gen-Alpha-typed-Integration scheme for standard constrained DAE
%
% - based only on constraint on position level
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - not derived from variational principle but simply evaluates RHS at
%   t_{n+1/2}
%
% Author: Philipp Kinon
% Date  : 09.12.2020

    methods
        
        function self = MP_std(this_simulation, this_system)
            self.DT    = this_simulation.DT;
            self.T_0   = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t     = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT    = size(self.t, 2) - 1;
            self.nVARS = 2*this_system.nDOF+this_system.mCONSTRAINTS;
            self.INDI_VELO = false;
            self.LM0   = zeros(this_system.mCONSTRAINTS,1);
            self.hasPARA = false;
            self.NAME  = 'MP-std';
        end
        
        function z0 = set_initial_condition(self,this_simulation,this_system)
            
            z0 = [this_simulation.Q_0', (this_system.MASS_MAT * this_simulation.V_0)', self.LM0'];
            
        end
        
        function [resi,tang] = compute_resi_tang(self,zn1,zn,this_system)
            
            %% Abbreviations
            M  = this_system.MASS_MAT;
            IM = M\eye(size(M));
            h  = self.DT;
            n  = this_system.nDOF;
            m  = this_system.mCONSTRAINTS;
            
            %% Unknows which will be iterated
            qn1       = zn1(1:n);
            pn1       = zn1(n+1:2*n);
            lambda_n1 = zn1(2*n+1:2*n+m);
            G_n1      = this_system.constraint_gradient(qn1);
            g_n1      = this_system.constraint(qn1);
            
            %% Known quantities from last time-step
            qn      = zn(1:n);
            pn      = zn(n+1:2*n);
            lambdan = zn(2*n+1:2*n+m);
            
            %% Alpha evaluated quantities
            q_n05      = 0.5*qn + 0.5*qn1;
            p_n05      = 0.5*pn + 0.5*pn1;
            lambda_n05 = 0.5*lambdan + 0.5*lambda_n1;
            DV_n05     = this_system.internal_potential_gradient(q_n05)+ this_system.external_potential_gradient(q_n05);
            D2V_n05    = this_system.internal_potential_hessian(q_n05) + this_system.external_potential_hessian(q_n05);
            G_n05      = this_system.constraint_gradient(q_n05);
            t_n05      = zeros(n);
            
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            for i = 1:m
                t_n05   = t_n05 + this_system.constraint_hessian(q_n05,m)*lambda_n05(m);
            end
            
            %% Residual vector 
            resi = [qn1 - qn - h*IM*p_n05                        ;
                    pn1 - pn + h*DV_n05   + h*G_n05'*lambda_n05  ;
                    g_n1                                         ];

            %% Tangent matrix
            tang = [eye(n)                                    -h*0.5*IM      zeros(n,m) ;
                    h*0.5*D2V_n05 + h*0.5*t_n05               eye(n)         0.5*h*G_n05'; 
                    G_n1                                      zeros(n,m)'    zeros(m)          ];
            
        end
        
    end

end