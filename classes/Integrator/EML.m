classdef EML < Integrator
    % Energy_Momentum-Integration scheme for standard constrained DAE
    %
    % - not derived from variational principle
    %
    % - uses Livens equations of motion
    %
    % - uses discrete gradient for ext. potential, internal potential and
    %   constraint gradients
    %
    % - G-equivarient versions for constraint and int. pot. gradients
    %
    % - takes account of non-constant mass-matrices
    %
    % Author: Philipp Kinon
    % Date  : 09.12.2022

    methods

        
        function self = EML(this_simulation, this_system)
        %% Constructor 
        
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF + this_system.mCONSTRAINTS;
            self.INDI_VELO = true;
            self.LM0 = zeros(this_system.mCONSTRAINTS, 1);
            self.hasPARA = false;
            self.NAME = 'EML';
            self.has_enhanced_constraint_force = [];
        end

        function z0 = set_initial_condition(self, this_simulation, this_system)
            p0 = (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0);
            z0 = [this_simulation.Q_0', p0', this_simulation.V_0', self.LM0'];
        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)
            % computes residual tangent
            %
            % :param zn1: input zn1
            % :param zn: input zn
            % :param this_system: input this_system
            % :returns: [ResidualVector, TangentMatrix]

            %% Abbreviations
            h = self.DT;
            nDOF = this_system.nDOF;
            mConstraints = this_system.mCONSTRAINTS;
            nPotInv = this_system.nPotentialInvariants;

            %% Unknows which will be iterated
            qn1 = zn1(1:nDOF);
            pn1 = zn1(nDOF+1:2*nDOF);
            vn1 = zn1(2*nDOF+1:3*nDOF);
            lambdan1 = zn1(3*nDOF+1:3*nDOF+mConstraints);
            Vext_n1 = this_system.external_potential(qn1);
            Mn1 = this_system.get_mass_matrix(qn1);
            g_n1 = this_system.constraint(qn1);

            %% Known quantities from last time-step
            qn = zn(1:nDOF);
            pn = zn(nDOF+1:2*nDOF);
            vn = zn(2*nDOF+1:3*nDOF);
            Vext_n = this_system.external_potential(qn);
            Mn = this_system.get_mass_matrix(qn);

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            v_n05 = 0.5 * (vn + vn1);
            DVext_n05 = this_system.external_potential_gradient(q_n05);
            D_1_T_n05 = this_system.kinetic_energy_gradient_from_velocity(q_n05, v_n05);
            
            % kinetic energy with mixed evaluations
            T_qn1vn  = 0.5 * vn'  * Mn1 * vn;
            T_qnvn   = 0.5 * vn'  * Mn  * vn;
            T_qn1vn1 = 0.5 * vn1' * Mn1 * vn1;
            T_qnvn1  = 0.5 * vn1' * Mn  * vn1;

            %% Discrete gradients

            % discrete gradient of the constraints
            DG_g = zeros(mConstraints, nDOF);
            g_invariants_difference_too_small = false;
        
            % for every constraint invariant individually
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
          

            % discrete gradient of internal potential based on invariants 
            DG_Vint = zeros(nDOF, 1); % for the internal potential
            V_invariants_difference_too_small = false;

            for i = 1:nPotInv % loop over all invariants
                %compute i-th invariants
                pi_n = this_system.potential_invariant(qn, i);
                pi_n1 = this_system.potential_invariant(qn1, i);
                % derivative of invariant w.r.t. q_n05
                DPiq_n05 = this_system.potential_invariant_gradient(q_n05, i);
                % evaluate internal potential depending on invariants
                Vs_n = this_system.potential_from_invariant(pi_n, i);
                Vs_n1 = this_system.potential_from_invariant(pi_n1, i);

                % if invariants at n and n1 are approx. equal use the midpoint
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
            
            % discrete gradients of kinetic energy and of external
            % potential energy
            D_1_T_qn05_vn = this_system.kinetic_energy_gradient_from_velocity(q_n05, vn);
            D_1_T_qn05_vn1 = this_system.kinetic_energy_gradient_from_velocity(q_n05, vn1);
            if abs((qn1-qn)'*(qn1-qn)) > 1e-9
                % discrete gradient of kinetic energy w.r.t position 
                DG_1_T_q_vn = D_1_T_qn05_vn + ((T_qn1vn - T_qnvn - D_1_T_qn05_vn'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_T_q_vn1 = D_1_T_qn05_vn1 + ((T_qn1vn1 - T_qnvn1 - D_1_T_qn05_vn1'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_T = 0.5*(DG_1_T_q_vn + DG_1_T_q_vn1);
                % discrete gradient of external potential energy
                DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn)) ) * (qn1-qn);
            else
                % use MP evaluation if qn1 is approx. qn
                DG_1_T = D_1_T_n05;
                DG_Vext = DVext_n05;
            end
            
            % discrete gradient of kinetic energy w.r.t velocity 
            DG_2_T = 0.5*(Mn + Mn1)*v_n05;

            %% Residual vector
            resi = [qn1 - qn - h * v_n05; 
                    pn1 - pn + h * DG_Vext + h * DG_Vint - h * DG_1_T + h * DG_g' * lambdan1; 
                    p_n05 - DG_2_T;
                    g_n1];

            %% Tangent matrix
            tang = [];

        end

    end

end
