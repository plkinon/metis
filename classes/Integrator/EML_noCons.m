classdef EML_noCons < Integrator

    %% Energy_Momentum-Integration scheme for standard constrained DAE
    %
    % - not derived from variational principle
    %
    % - uses Livens equations of motion
    %
    % - uses standard gradient for ext. potential and discrete gradient for
    %   internal potential
    %
    % - takes account of non-constant mass-matrices
    %
    % Author: Philipp Kinon
    % Date  : 09.12.2022

    methods

        function self = EML_noCons(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF;
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'EMG-noCons';
            self.has_enhanced_constraint_force = [];
end

        function z0 = set_initial_condition(~, this_simulation, this_system)

            z0 = [this_simulation.Q_0', (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0)', this_simulation.V_0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)

            %% Abbreviations
            h = self.DT;
            n = this_system.nDOF;
            p = this_system.nPotentialInvariants;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);
            Vext_n1 = this_system.external_potential(qn1);
            Mn1 = this_system.get_mass_matrix(qn1);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            vn = zn(2*n+1:3*n);
            Vext_n = this_system.external_potential(qn);
            Mn = this_system.get_mass_matrix(qn);

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            v_n05 = 0.5 * (vn + vn1);
            DVext_n05 = this_system.external_potential_gradient(q_n05);
            D2Vext_n05 = this_system.external_potential_hessian(q_n05);
            D2Vint_n05 = this_system.internal_potential_hessian(q_n05);
            D_1_T_n05 = this_system.kinetic_energy_gradient_from_velocity(q_n05, v_n05);

            %% Discrete gradients
            T_qn1vn  = 0.5 * vn'  * Mn1 * vn;
            T_qnvn   = 0.5 * vn'  * Mn  * vn;
            T_qn1vn1 = 0.5 * vn1' * Mn1 * vn1;
            T_qnvn1  = 0.5 * vn1' * Mn  * vn1;
            % for the internal potential
            DG_Vint = zeros(n, 1); % for the internal potential
            DG_Vext = zeros(n, 1); % for the external potential
            DG_1_T = zeros(n, 1);
            D_1_T_qn05_vn = this_system.kinetic_energy_gradient_from_velocity(q_n05, vn);
            D_1_T_qn05_vn1 = this_system.kinetic_energy_gradient_from_velocity(q_n05, vn1);
            K21_DG_V = zeros(n, n);
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

                %for the tangent matrix
                D2PiDq2 = this_system.potential_invariant_hessian(q_n05, i);
                DPiDq_n1 = this_system.potential_invariant_gradient(qn1, i);
                DVsDpi_n1 = this_system.potential_gradient_from_invariant(pi_n1, i);

                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(pi_n1-pi_n) > 1e-09
                    % discrete gradient
                    DG_Vint = DG_Vint + (Vs_n1 - Vs_n) / (pi_n1 - pi_n) * DPiq_n05;
                    K21_DG_V = K21_DG_V + DPiq_n05 * (DVsDpi_n1 * DPiDq_n1 * 1 / (pi_n1 - pi_n) - (Vs_n1 - Vs_n) / (pi_n1 - pi_n)^2 * DPiDq_n1)' + (Vs_n1 - Vs_n) / (pi_n1 - pi_n) * 1 / 2 * D2PiDq2;
                else
                    V_invariants_difference_too_small = true;
                    break
                end

            end

            if V_invariants_difference_too_small
                % else use MP evaluation of gradient
                DG_Vint = this_system.internal_potential_gradient(q_n05);
                K21_DG_V = 1 / 2 * D2Vint_n05;
            end

            if abs((qn1-qn)'*(qn1-qn)) > 1e-9
                DG_1_T_q_vn = D_1_T_qn05_vn + ((T_qn1vn - T_qnvn - D_1_T_qn05_vn'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_T_q_vn1 = D_1_T_qn05_vn1 + ((T_qn1vn1 - T_qnvn1 - D_1_T_qn05_vn1'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_T = 0.5*(DG_1_T_q_vn + DG_1_T_q_vn1);
                DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn)) ) * (qn1-qn);
            else
                DG_1_T = D_1_T_n05;
                DG_Vext = DVext_n05;
            end

            DG_2_T = 0.5*(Mn + Mn1)*v_n05;

            %% Residual vector
            resi = [qn1 - qn - h * v_n05; 
                    pn1 - pn + h * DG_Vext + h * DG_Vint - h * DG_1_T;
                    p_n05 - DG_2_T];

            %% Tangent matrix
            %tang = [eye(n), -h * 0.5 * IM; 
            %        h * 0.5 * D2Vext_n05 + h * K21_DG_V, eye(n)];
            tang = [];
        end

    end

end