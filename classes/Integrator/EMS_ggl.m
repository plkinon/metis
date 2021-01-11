classdef EMS_ggl < Integrator
%% Energy_Momentum-Integration scheme for constrained DAE with GGL stabilisation
%
% - based only on constraints on position and velocity level
%
% - independent momenta variables (Hamilton Potryagin approach)
%
% - not derived from variational principle
% 
% - basic approach from Gonzales 1999 here also applied to constraint on
%   velocity level
%
% - uses standard gradient for ext. potential and discrete gradient for
%   internal potential and constraints
%
% Author: Philipp Kinon
% Date  : 20.12.2020

    methods
        function self = initialise(self,~,this_problem)
            self.NAME  = 'EMS-GGL';
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
            p  = this_problem.nPotentialInvariants;
            
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
            
            %% MP evaluated quantities
            q_n05      = 0.5*(qn + qn1);
            p_n05      = 0.5*(pn + pn1);
            DVext_n05  = this_problem.external_potential_gradient(q_n05);
            D2Vext_n05 = this_problem.external_potential_hessian(q_n05);
            D2Vint_n05 = this_problem.internal_potential_hessian(q_n05);
            t_n05l      = zeros(n);
            t_n05g      = zeros(n);
            
            %% Discrete gradients
            % for the internal potential
            DG_Vint = zeros(n,1);
            K21_DG_V = zeros(n,n);
            
            % for every invariant individually
            V_invariants_difference_too_small = false;
            for i = 1:p
                %compute i-th invariants 
                pi_n     = this_problem.potential_invariant(qn,i);
                pi_n1    = this_problem.potential_invariant(qn1,i);
                % evaluate internal potential depending on invariants
                Vs_n     = this_problem.potential_from_invariant(pi_n,i);
                Vs_n1    = this_problem.potential_from_invariant(pi_n1,i);
                % derivative of invariant w.r.t. q_n05
                DPiq_n05 = this_problem.potential_invariant_gradient(q_n05,i);
                
                %for the tangent matrix
                D2PiDq2   = this_problem.potential_invariant_hessian(q_n05,i);
                DPiDq_n1  = this_problem.potential_invariant_gradient(qn1,i);
                DVsDpi_n1 = this_problem.potential_gradient_from_invariant(pi_n1,i);
                
                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(pi_n1 - pi_n) > 1e-09
                    
                    % discrete gradient
                    DG_Vint  = DG_Vint + (Vs_n1-Vs_n)/(pi_n1-pi_n)*DPiq_n05;
                    K21_DG_V    = K21_DG_V + DPiq_n05*(DVsDpi_n1*DPiDq_n1*1/(pi_n1-pi_n) - (Vs_n1-Vs_n)/(pi_n1-pi_n)^2*DPiDq_n1)' + (Vs_n1-Vs_n)/(pi_n1-pi_n)*1/2*D2PiDq2;
                    
                else
                    V_invariants_difference_too_small = true;
                    break
                end
                
            end
            
            if V_invariants_difference_too_small
                 % else use MP evaluation of gradient
                    DG_Vint = this_problem.internal_potential_gradient(q_n05);
                    K21_DG_V  = K21_DG_V  + 1/2*D2Vint_n05;
            end
            
            %% for the gradients of the position constraints
            DG_g    = zeros(m,n);
            K21_DG_g = zeros(n,n);
            
            % for every invariant individually
            for j = 1:this_problem.nConstraintInvariants
                %compute i-th invariants 
                zeta_n  = this_problem.constraint_invariant(qn,j);
                zeta_n1 = this_problem.constraint_invariant(qn1,j);
                % evaluate constraints depending on invariants
                gs_n = this_problem.constraint_from_invariant(zeta_n,j);
                gs_n1 = this_problem.constraint_from_invariant(zeta_n1,j);
                % derivative of invariant w.r.t. q_n05
                DzetaDq_n05 = this_problem.constraint_invariant_gradient(q_n05,j);
                
                % tangent matrix terms
                D2zetaDq2   = this_problem.constraint_invariant_hessian(qn1,j);
                DgsDzeta_n1 = this_problem.constraint_gradient_from_invariant(qn1,j);
                DzetaDq_n1  = this_problem.constraint_invariant_gradient(qn1,j);
                
                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(zeta_n1 - zeta_n) > 1e-9
                    
                    % discrete gradient for poisition constraint
                    DG_g(j,:) = (gs_n1 - gs_n)/(zeta_n1 - zeta_n)*DzetaDq_n05';
                    
                    K21_DG_g    = K21_DG_g + lambda_n1(j)*(DzetaDq_n05'*(DgsDzeta_n1*DzetaDq_n1*1/(zeta_n1-zeta_n) - (gs_n1 - gs_n)/(zeta_n1 - zeta_n)^2*DzetaDq_n1) + (gs_n1 - gs_n)/(zeta_n1 - zeta_n)*1/2*D2zetaDq2);
                    g_invariants_difference_too_small = false;
                    
                else
                    g_invariants_difference_too_small = true;
                    break
                end
            end
            
            if g_invariants_difference_too_small
                % else use MP evaluation of gradient
                G_n05 = this_problem.constraint_gradient(q_n05);
                DG_g  = G_n05;
                for j=1:this_problem.nConstraintInvariants
                    D2g_Dq_n05 = this_problem.constraint_hessian(q_n05,j);
                    K21_DG_g  = K21_DG_g  + 1/2*lambda_n1(j)*D2g_Dq_n05;
                end
            end
            
            %% for the gradients of the velocity constraints
            DG_gv_q = zeros(m,n);
            DG_gv_p = zeros(m,n);
            K11_DG_gv = zeros(n,n);
            K12_DG_gv = zeros(n,n);
            K21_DG_gv = zeros(n,n);
            K22_DG_gv = zeros(n,n);
            
            % for every invariant individually
            for k = 1:this_problem.nVconstraintInvariants
                
                %compute k-th invariants 
                pi_n = this_problem.vConstraint_invariant(qn,pn,k);
                pi_n1 = this_problem.vConstraint_invariant(qn1,pn1,k);
                % evaluate constraints depending on invariants
                gv_n1 = this_problem.Vconstraint_from_invariant(pi_n1,k);
                gv_n  = this_problem.Vconstraint_from_invariant(pi_n,k);
                % derivative of invariant w.r.t. q_n05
                DpiDqn05 = this_problem.vConstraint_invariant_gradient_q(q_n05,p_n05,k);
                DpiDpn05 = this_problem.vConstraint_invariant_gradient_p(q_n05,p_n05,k);
                
                % tangent matrix terms
                DgvDpi_n1 = this_problem.Vconstraint_gradient_from_invariant(pi_n1,k);
                DpiDq_n1 = this_problem.vConstraint_invariant_gradient_q(qn1,pn1,k);
                DpiDp_n1 = this_problem.vConstraint_invariant_gradient_p(qn1,pn1,k);
                D2piDqDp = this_problem.vConstraint_invariant_hessian_qp(q_n05,p_n05,k);
                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(pi_n1 - pi_n) > 1e-09
                    
                    % discrete gradient for velocity constraint
                    DG_gv_q(k,:) = (gv_n1 -gv_n)/(pi_n1-pi_n)*DpiDqn05;
                    DG_gv_p(k,:) = (gv_n1 -gv_n)/(pi_n1-pi_n)*DpiDpn05;
                    
                    K11_DG_gv = K11_DG_gv + gamma_n1(k)*(DpiDpn05*(DgvDpi_n1*DpiDq_n1*1/(pi_n1-pi_n) - (gv_n1-gv_n)/(pi_n1-pi_n)^2*DpiDq_n1)' + (gv_n1-gv_n)/(pi_n1-pi_n)*1/2*D2piDqDp);
                    K12_DG_gv = K12_DG_gv + gamma_n1(k)*(DpiDpn05*(DgvDpi_n1*DpiDp_n1*1/(pi_n1-pi_n) - (gv_n1-gv_n)/(pi_n1-pi_n)^2*DpiDp_n1)' + (gv_n1-gv_n)/(pi_n1-pi_n)*1/2*0);
                    K21_DG_gv = K21_DG_gv + gamma_n1(k)*(DpiDqn05*(DgvDpi_n1*DpiDq_n1*1/(pi_n1-pi_n) - (gv_n1-gv_n)/(pi_n1-pi_n)^2*DpiDq_n1)' + (gv_n1-gv_n)/(pi_n1-pi_n)*1/2*0);
                    K22_DG_gv = K22_DG_gv + gamma_n1(k)*(DpiDqn05*(DgvDpi_n1*DpiDp_n1*1/(pi_n1-pi_n) - (gv_n1-gv_n)/(pi_n1-pi_n)^2*DpiDp_n1)' + (gv_n1-gv_n)/(pi_n1-pi_n)*1/2*D2piDqDp);
                    gv_invariants_difference_too_small = false;
                else
                    gv_invariants_difference_too_small = true;
                    break
                    % else use MP evaluation of gradient
                    
                end
            end
            
            if gv_invariants_difference_too_small
                % else use MP evaluation of gradient
                for k = 1:this_problem.nVconstraintInvariants
                    D2g_qq_tmp   = this_problem.constraint_hessian(q_n05,k);
                    G_n05        = this_problem.constraint_gradient(q_n05);
                    DG_gv_q(k,:) = D2g_qq_tmp*IM*p_n05;
                    DG_gv_p(k,:) = G_n05(k,:)*IM;
                    K11_DG_gv    = K11_DG_gv + gamma_n1(k)*1/2*this_problem.constraint_hessian(q_n05,k)*IM;
                    K22_DG_gv    = K22_DG_gv + gamma_n1(k)*1/2*D2g_qq_tmp*IM;
                end
                K12_DG_gv = zeros(n);
                K21_DG_gv = zeros(n);
            end
            
            %% Terms with Hessian of constraints
            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            for i = 1:m
                t_n05l   = t_n05l + this_problem.constraint_hessian(q_n05,m)*lambda_n1(m);
                t_n05g   = t_n05g + this_problem.constraint_hessian(q_n05,m)*gamma_n1(m);
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
            resi = [qn1 - qn - h*IM*p_n05    - h*DG_gv_p'*gamma_n1                                    ;
                    pn1 - pn + h*DVext_n05   + h*DG_Vint + h*DG_g'*lambda_n1 + h*DG_gv_q'*gamma_n1 ;
                    g_n1                                                                           ;
                    G_n1*IM*pn1                                                                   ];

            %% Tangent matrix
            %tang = [];
            tang = [eye(n)-h*K11_DG_gv                                          -h*0.5*IM+h*K12_DG_gv         zeros(n,m)    -h*DG_gv_p';
                    h*0.5*D2Vext_n05 + h*K21_DG_V + h*K21_DG_g + h*K21_DG_gv    eye(n) + h*K22_DG_gv          h*DG_g'       h*DG_gv_q'       ; 
                    G_n1                                                        zeros(n,m)'                 zeros(m)      zeros(m)         ;
                    T_n1                                                        G_n1*IM                     zeros(m)      zeros(m)      ];
        end
        
    end

end