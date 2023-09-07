classdef EMS_std < Integrator

    %% Energy_Momentum-Integration scheme for standard constrained DAE
    %
    % - based only on constraint on position level
    %
    % - independent momenta variables (Hamilton Potryagin approach)
    %
    % - not derived from variational principle
    %
    % - taken from Gonzales 1999
    %
    % - uses standard gradient for ext. potential and discrete gradient for
    %   internal potential and constraint
    %
    % Author: Philipp Kinon
    % Date  : 17.12.2020

    methods

        function self = EMS_std(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 2 * this_system.nDOF + 1 * this_system.mCONSTRAINTS;
            self.INDI_VELO = false;
            self.LM0 = zeros(this_system.mCONSTRAINTS, 1);
            self.hasPARA = false;
            self.NAME = 'EMS-std';
            self.has_enhanced_constraint_force = false;
        end

        function z0 = set_initial_condition(self, this_simulation, this_system)
            q0 = this_simulation.Q_0;
            M0 = this_system.get_mass_matrix(q0);
            v0 = this_simulation.V_0;
            z0 = [q0', (M0 * v0)', self.LM0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)

            %% Abbreviations
            h = self.DT;
            nDOF = this_system.nDOF;
            m = this_system.mCONSTRAINTS;
            p = this_system.nPotentialInvariants;
            nKinInv = this_system.nKineticInvariants;

            %% Unknows which will be iterated
            qn1 = zn1(1:nDOF);
            pn1 = zn1(nDOF+1:2*nDOF);
            lambda_n1 = zn1(2*nDOF+1:2*nDOF+m);
            g_n1 = this_system.constraint(qn1);
            Mn1 = this_system.get_mass_matrix(qn1);
            IMn1 = eye(size(Mn1)) / Mn1;

            %% Known quantities from last time-step
            qn = zn(1:nDOF);
            pn = zn(nDOF+1:2*nDOF);
            Mn = this_system.get_mass_matrix(qn);
            IMn = eye(size(Mn)) / Mn;

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            Mn05 = this_system.get_mass_matrix(q_n05);
            IMn05 = eye(size(Mn05)) / Mn05;

            DVext_n05 = this_system.external_potential_gradient(q_n05);
            D_1_T_n05 = this_system.kinetic_energy_gradient_from_momentum(q_n05, p_n05);

            % kinetic energy with mixed evaluations
            T_qn1pn  = 0.5 * pn'  * IMn1 * pn;
            T_qnpn   = 0.5 * pn'  * IMn  * pn;
            T_qn1pn1 = 0.5 * pn1' * IMn1 * pn1;
            T_qnpn1  = 0.5 * pn1' * IMn  * pn1;

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
            end
            
            % kinetic energy discrete gradient
            DG_T_q = zeros(nDOF, 1); % for the kinetic energy
            DG_T_p = zeros(nDOF, 1); % for the kinetic energy
            T_invariants_difference_too_small = false;

            if  ~any(nKinInv)
                   % discrete gradients of kinetic energy and of external
                    % potential energy
                    D_1_T_qn05_pn = this_system.kinetic_energy_gradient_from_momentum(q_n05, pn);
                    D_1_T_qn05_pn1 = this_system.kinetic_energy_gradient_from_momentum(q_n05, pn1);
                    if abs((qn1-qn)'*(qn1-qn)) > 1e-9
                        % discrete gradient of kinetic energy w.r.t position 
                        DG_1_T_q_pn = D_1_T_qn05_pn + ((T_qn1pn - T_qnpn - D_1_T_qn05_pn'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                        DG_1_T_q_pn1 = D_1_T_qn05_pn1 + ((T_qn1pn1 - T_qnpn1 - D_1_T_qn05_pn1'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                        DG_T_q = 0.5*(DG_1_T_q_pn + DG_1_T_q_pn1);
                    else
                        % use MP evaluation if qn1 is approx. qn
                        DG_T_q = D_1_T_n05;
                    end
                    
                    % discrete gradient of kinetic energy w.r.t velocity 
                    DG_T_p = 0.5*(IMn + IMn1)*p_n05;
            else


                for k = 1:nKinInv % loop over all quadratic invariants
                    %compute i-th invariants
                    omega_n = this_system.kinetic_energy_invariant(qn, pn, k);
                    omega_n1 = this_system.kinetic_energy_invariant(qn1, pn1, k);
                    omega_n05 = 1/2*(omega_n+omega_n1);
                    % derivative of invariant w.r.t. q_n05
                    Domegaq_n05 = this_system.kinetic_energy_invariant_gradient_q(q_n05, p_n05, k);
                    % derivative of invariant w.r.t. v_n05
                    Domegap_n05 = this_system.kinetic_energy_invariant_gradient_p(q_n05, p_n05, k);
                    % evaluate internal potential depending on invariants
                    Ts_n = this_system.kinetic_energy_from_invariant(omega_n, k);
                    Ts_n1 = this_system.kinetic_energy_from_invariant(omega_n1, k);
    
                    % if invariants at n and n1 are approx. equal use the midpoint
                    % evaluated gradient instead
                    if abs((omega_n1-omega_n)'*(omega_n1-omega_n)) > 1e-09
                        % discrete gradient
                        DT_omega_n05 = this_system.kinetic_energy_gradient_from_invariant(omega_n05,k);
                        DT_omega = DT_omega_n05 + ((Ts_n1 - Ts_n - DT_omega_n05'*(omega_n1 -omega_n)) / ((omega_n1-omega_n)'*(omega_n1-omega_n))) * (omega_n1-omega_n); 
                        % if the kinetic energy is quadratic in this invariant, the second term vanishes
                        DG_T_q = DG_T_q + Domegaq_n05' * DT_omega;
                        DG_T_p = DG_T_p + Domegap_n05' * DT_omega;

                    else
                        T_invariants_difference_too_small = true;
                        break
                    end
    
                end
    
                if T_invariants_difference_too_small
                    % else use MP evaluation of gradient
                    DG_T_q = this_system.kinetic_energy_gradient_from_momentum(q_n05, p_n05);
                    DG_T_p = IMn05*p_n05;
                end
            

            end

            % for the gradients of the constraints
            DG_g = zeros(m, nDOF);
            K21_DG_g = zeros(nDOF, nDOF);
            g_invariants_difference_too_small = false;

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

                % tangent matrix terms
                D2zetaDq2 = this_system.constraint_invariant_hessian(qn1, j);
                DgsDzeta_n1 = this_system.constraint_gradient_from_invariant(qn1, j);
                DzetaDq_n1 = this_system.constraint_invariant_gradient(qn1, j);

                % if invariants at n and n1 are equal use the midpoint
                % evaluated gradient instead
                if abs(zeta_n1-zeta_n) > 1e-9
                    % discrete gradient
                    DG_g(j, :) = (gs_n1 - gs_n) / (zeta_n1 - zeta_n) * DzetaDq_n05';
                    K21_DG_g = K21_DG_g + lambda_n1(j) * (DzetaDq_n05' * (DgsDzeta_n1 * DzetaDq_n1 * 1 / (zeta_n1 - zeta_n) - (gs_n1 - gs_n) / (zeta_n1 - zeta_n)^2 * DzetaDq_n1) + (gs_n1 - gs_n) / (zeta_n1 - zeta_n) * 1 / 2 * D2zetaDq2);
                else
                    g_invariants_difference_too_small = true;
                    break
                end

            end

            if g_invariants_difference_too_small
                % else use MP evaluation of gradient
                G_n05 = this_system.constraint_gradient(q_n05);
                DG_g = G_n05;
                for j = 1:this_system.nConstraintInvariants
                    D2g_Dq_n05 = this_system.constraint_hessian(q_n05, j);
                    K21_DG_g = K21_DG_g + 1 / 2 * lambda_n1(j) * D2g_Dq_n05;
                end
            end

            %% Residual vector
% for constant mass matrix:
%             resi = [qn1 - qn - h * IM * p_n05; 
%                     pn1 - pn + h * DVext_n05 + h * DG_Vint + h * DG_g' * lambda_n1; 
%                     g_n1];
            resi = [qn1 - qn - h * DG_T_p; 
                    pn1 - pn + h * DG_T_q + h * DVext_n05 + h * DG_Vint + h * DG_g' * lambda_n1; 
                    g_n1];

            %% Tangent matrix
            tang = [];

        end

    end

end