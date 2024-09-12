classdef PHDAE_DG < Integrator
    % Discrete gradient Integration scheme for port-Hamiltonian
    % semiexplicit DAEs
    %
    % - based only on constraint on velocity level
    %
    % - not derived from variational principle
    %
    % - assumes constant mass matrix
    
    methods
        
        function self = PHDAE_DG(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 2 * this_system.nDOF + 1 * this_system.mCONSTRAINTS;
            self.INDI_VELO = true;
            self.LM0 = zeros(this_system.mCONSTRAINTS, 1);
            self.hasPARA = false;
            self.NAME = 'PHDAE-DG';
            self.has_enhanced_constraint_force = false;
        end
        
        function z0 = set_initial_condition(self, this_simulation, this_system)
            q0 = this_simulation.Q_0;
            M0 = this_system.get_mass_matrix(q0);
            v0 = this_simulation.V_0;
            z0 = [q0', v0', self.LM0'];
            
        end
        
        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system, time_n)
            % Computes residual vector & tangent matrix
            %
            % :param zn1: state vector for next time step
            % :param zn: state vector at current time step
            % :param this_system: System object
            % :returns: [ResidualVector, TangentMatrix] for the Newton's method to update zn1
            
            %% Abbreviations
            h = self.DT;
            nDOF = this_system.nDOF;
            m = this_system.mCONSTRAINTS;
            p = this_system.nPotentialInvariants;
            
            %% Unknows which will be iterated
            qn1 = zn1(1:nDOF);
            vn1 = zn1(nDOF+1:2*nDOF);
            lambda_n1 = zn1(2*nDOF+1:2*nDOF+m);
            Vext_n1 = this_system.external_potential(qn1);
            g_n1 = this_system.constraint(qn1);
            M = this_system.get_mass_matrix(qn1);
            
            %% Known quantities from last time-step
            qn = zn(1:nDOF);
            vn = zn(nDOF+1:2*nDOF);
            Vext_n = this_system.external_potential(qn);
            
            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            v_n05 = 0.5 * (vn + vn1);
            DVext_n05 = this_system.external_potential_gradient(q_n05);
            
            
            % for the gradients of the constraints
            DG_g = zeros(m, nDOF);
            %K21_DG_g = zeros(nDOF, nDOF);
            g_invariants_difference_too_small = false;

            if ismethod(this_system,'get_dissipation_matrix')
                Dn05 = this_system.get_dissipation_matrix(q_n05);
            end
            
            % for every invariant individually
            for j = 1:this_system.nConstraintInvariants
                %compute i-th invariants
                zeta_n = this_system.constraint_invariant(qn, j);
                zeta_n1 = this_system.constraint_invariant(qn1, j);
                
                % evaluate constraints depending on invariants
                gs_n = this_system.constraint_from_invariant(zeta_n, j);
                gs_n1 = this_system.constraint_from_invariant(zeta_n1, j);
                % derivative of invariant w.r.t. q_n05
                DzetaDq_n05 = this_system.constraint_invariant_gradient(q_n05, j);
                
                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(zeta_n1-zeta_n) > 1e-9
                    % discrete gradient
                    DG_g(j, :) = (gs_n1 - gs_n) / (zeta_n1 - zeta_n) * DzetaDq_n05';
                else
                    g_invariants_difference_too_small = true;
                    break
                end
                
            end
            
            if g_invariants_difference_too_small
                % else use MP evaluation of gradient
                G_n05 = this_system.constraint_gradient(q_n05);
                DG_g = G_n05;
            end
            
            %% Discrete gradients
            % for the internal potential
            DG_Vint = zeros(nDOF, 1);
            V_invariants_difference_too_small = false;
            
            % for every invariant individually
            for i = 1:p
                %compute i-th invariants
                pi_n = this_system.potential_invariant(qn, i);
                pi_n1 = this_system.potential_invariant(qn1, i);
                % derivative of invariant w.r.t. q_n05
                DPiq_n05 = this_system.potential_invariant_gradient(q_n05, i);
                % evaluate internal potential depending on invariants
                Vs_n = this_system.potential_from_invariant(pi_n, i);
                Vs_n1 = this_system.potential_from_invariant(pi_n1, i);
                
                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(pi_n1-pi_n) > 1e-09
                    % discrete gradient
                    DG_Vint = DG_Vint + (Vs_n1 - Vs_n) / (pi_n1 - pi_n) * DPiq_n05;
                else
                    V_invariants_difference_too_small = true;
                    break
                end
                
            end
            
            if V_invariants_difference_too_small
                % else use MP evaluation of gradient
                DG_Vint = this_system.internal_potential_gradient(q_n05);
            end
            
            %% Potential energy gradient
            if abs((qn1-qn)'*(qn1-qn)) > 1e-9
                % discrete gradient of external potential energy
                DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn)) ) * (qn1-qn);
            else
                DG_Vext = DVext_n05;
            end

            
            %% Residual vector
            resi = [qn1 - qn - h * v_n05;
                    M*(vn1 - vn) + h * DG_Vext + h * DG_Vint + h * Dn05 * v_n05 + h * DG_g' * lambda_n1;
                    DG_g*v_n05];
            
            %% Tangent matrix
            tang = [];
            
        end
        
    end
    
end