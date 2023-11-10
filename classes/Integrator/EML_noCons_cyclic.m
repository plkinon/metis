classdef EML_noCons_cyclic < Integrator
    % Energy_Momentum-Integration scheme for ODE with cyclic coordinate
    %
    % - not derived from variational principle
    %
    % - uses Livens equations of motion
    %
    % - uses discrete gradient for ext. potential and internal potential
    %
    % - takes account of non-constant mass-matrices
    %
    % - splits up for cyclic coordinates! import for conservation of
    %   corresponding conjugate momenta


    methods

        function self = EML_noCons_cyclic(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF;
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'EML-noCons';
            self.has_enhanced_constraint_force = [];
end

        function z0 = set_initial_condition(~, this_simulation, this_system)

            z0 = [this_simulation.Q_0', (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0)', this_simulation.V_0'];

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
            xn1 = qn1(~this_system.isCyclicCoordinate);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);
            Vext_n1 = this_system.external_potential_cyclic(xn1);
            Mn1 = this_system.get_mass_matrix_cyclic(xn1);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            vn = zn(2*n+1:3*n);
            xn = qn(~this_system.isCyclicCoordinate);
            Vext_n = this_system.external_potential_cyclic(xn);
            Mn = this_system.get_mass_matrix_cyclic(xn);

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            v_n05 = 0.5 * (vn + vn1);
            x_n05 = 0.5 * (xn + xn1);
            DVext_n05 = this_system.external_potential_gradient_cyclic(x_n05);
            D_1_T_n05 = this_system.kinetic_energy_gradient_from_velocity_cyclic(x_n05, v_n05);

            %% Discrete gradients
            T_xn1vn  = 0.5 * vn'  * Mn1 * vn;
            T_xnvn   = 0.5 * vn'  * Mn  * vn;
            T_xn1vn1 = 0.5 * vn1' * Mn1 * vn1;
            T_xnvn1  = 0.5 * vn1' * Mn  * vn1;
            
            % for the internal potential
            DG_Vint = zeros(n, 1); % for the internal potential
            D_1_T_xn05_vn = this_system.kinetic_energy_gradient_from_velocity_cyclic(x_n05, vn);
            D_1_T_xn05_vn1 = this_system.kinetic_energy_gradient_from_velocity_cyclic(x_n05, vn1);
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
                if abs(pi_n1-pi_n) > 1e-07
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

            if abs((xn1-xn)'*(xn1-xn)) > 1e-9
                DG_1_T_x_vn = D_1_T_xn05_vn + ((T_xn1vn - T_xnvn - D_1_T_xn05_vn'*(xn1 -xn)) / ((xn1-xn)'*(xn1-xn))) * (xn1-xn); 
                DG_1_T_x_vn1 = D_1_T_xn05_vn1 + ((T_xn1vn1 - T_xnvn1 - D_1_T_xn05_vn1'*(xn1 -xn)) / ((xn1-xn)'*(xn1-xn))) * (xn1-xn);
                DG_1_T = 0.5*(DG_1_T_x_vn + DG_1_T_x_vn1);

                DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(xn1 -xn)) / ((xn1-xn)'*(xn1-xn))) * (xn1-xn);
            else
                DG_1_T = D_1_T_n05;
                DG_Vext = DVext_n05;
            end

            DG_2_T = 0.5*(Mn + Mn1)*v_n05;
            DG_1_T = [DG_1_T;0];
            DG_Vext = [DG_Vext;0];

            

            %% Residual vector
            resi = [qn1 - qn - h * v_n05; 
                    pn1 - pn + h * DG_Vext + h * DG_Vint - h * DG_1_T;
                    p_n05 - DG_2_T];

            %% Tangent matrix
            tang = [];
        end

    end

end