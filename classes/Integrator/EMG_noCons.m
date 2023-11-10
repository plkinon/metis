classdef EMG_noCons < Integrator
    % Energy_Momentum-Integration scheme for ODE
    %
    % - not derived from variational principle
    %
    % - taken from Gonzales 1996
    %
    % - uses standard midpoint gradient for ext. potential and discrete gradient for
    %   internal potential
    %
    % - takes account of non-constant mass-matrices


    methods

        function self = EMG_noCons(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 2 * this_system.nDOF;
            self.INDI_VELO = false;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'EMG-noCons';
            self.has_enhanced_constraint_force = [];
end

        function z0 = set_initial_condition(~, this_simulation, this_system)

            z0 = [this_simulation.Q_0', (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0)'];

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
            n = this_system.nDOF;
            p = this_system.nPotentialInvariants;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            Vext_n1 = this_system.external_potential(qn1);
            Mn1 = this_system.get_mass_matrix(qn1);
            IMn1 = Mn1 \ eye(size(Mn1));

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            Vext_n = this_system.external_potential(qn);
            Mn = this_system.get_mass_matrix(qn);
            IMn = Mn \ eye(size(Mn));

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            DVext_n05 = this_system.external_potential_gradient(q_n05);
            D2Vext_n05 = this_system.external_potential_hessian(q_n05);
            D2Vint_n05 = this_system.internal_potential_hessian(q_n05);
            D_1_T_n05 = this_system.kinetic_energy_gradient_from_momentum(q_n05, p_n05);
            D_1_T_qn05_pn = this_system.kinetic_energy_gradient_from_momentum(q_n05, pn);
            D_1_T_qn05_pn1 = this_system.kinetic_energy_gradient_from_momentum(q_n05, pn1);
            %% Discrete gradients
            T_qn1pn  = 0.5 * pn'  * IMn1 * pn;
            T_qnpn   = 0.5 * pn'  * IMn  * pn;
            T_qn1pn1 = 0.5 * pn1' * IMn1 * pn1;
            T_qnpn1  = 0.5 * pn1' * IMn  * pn1;
            DG_Vint = zeros(n, 1); % for the internal potential
            DG_Vext = zeros(n, 1); % for the external potential
            DG_1_T = zeros(n, 1); % for the kinetic energy
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

            if abs((qn1-qn)'*(qn1-qn)) > 1e-09
                DG_1_T_q_pn = D_1_T_qn05_pn + ((T_qn1pn - T_qnpn - D_1_T_qn05_pn'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_T_q_pn1 = D_1_T_qn05_pn1 + ((T_qn1pn1 - T_qnpn1 - D_1_T_qn05_pn1'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_T = 0.5*(DG_1_T_q_pn + DG_1_T_q_pn1);
                DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn)) ) * (qn1-qn);
            else
                DG_1_T = D_1_T_n05;
                DG_Vext = DVext_n05;
            end

            DG_2_T = 0.5*(IMn + IMn1)*p_n05;

            %% Residual vector
            resi = [qn1 - qn - h * DG_2_T; 
                    pn1 - pn + h * DG_Vext + h * DG_Vint + h * DG_1_T];

            %% Tangent matrix
            %tang = [eye(n), -h * 0.5 * IM; 
            %        h * 0.5 * D2Vext_n05 + h * K21_DG_V, eye(n)];
            tang = [];
        end

    end

end