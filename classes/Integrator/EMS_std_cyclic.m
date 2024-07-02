classdef EMS_std_cyclic < Integrator

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

        function self = EMS_std_cyclic(this_simulation, this_system)
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
            n = this_system.nDOF;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            xn1 = qn1(~this_system.isCyclicCoordinate);
            pn1 = zn1(n+1:2*n);
            Mn1 = this_system.get_mass_matrix_cyclic(xn1);
            IMn1 = eye(size(Mn1)) / Mn1;
            Vext_n1 = this_system.external_potential_cyclic(xn1);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            xn = qn(~this_system.isCyclicCoordinate);
            Vext_n = this_system.external_potential_cyclic(xn);
            Mn = this_system.get_mass_matrix_cyclic(xn);
            IMn = eye(size(Mn)) / Mn;

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            x_n05 = 0.5 * (xn + xn1);
            DVext_n05 = this_system.external_potential_gradient_cyclic(x_n05);
            Mn05 = this_system.get_mass_matrix(q_n05);
            IMn05 = eye(size(Mn05)) / Mn05;

            D_1_T_n05 = this_system.kinetic_energy_gradient_from_momentum_cyclic(x_n05, p_n05);

            % kinetic energy with mixed evaluations
            T_xn1pn  = 0.5 * pn'  * IMn1 * pn;
            T_xnpn   = 0.5 * pn'  * IMn  * pn;
            T_xn1pn1 = 0.5 * pn1' * IMn1 * pn1;
            T_xnpn1  = 0.5 * pn1' * IMn  * pn1;

            %% Discrete gradients

            % discrete gradients of kinetic energy and of external
            % potential energy
            D_1_T_xn05_pn = this_system.kinetic_energy_gradient_from_momentum_cyclic(x_n05, pn);
            D_1_T_xn05_pn1 = this_system.kinetic_energy_gradient_from_momentum_cyclic(x_n05, pn1);
            if abs((xn1-xn)'*(xn1-xn)) > 1e-9
                % discrete gradient of kinetic energy w.r.t position 
                DG_1_T_x_pn = D_1_T_xn05_pn + ((T_xn1pn - T_xnpn - D_1_T_xn05_pn'*(xn1 -xn)) / ((xn1-xn)'*(xn1-xn))) * (xn1-xn); 
                DG_1_T_x_pn1 = D_1_T_xn05_pn1 + ((T_xn1pn1 - T_xnpn1 - D_1_T_xn05_pn1'*(xn1 -xn)) / ((xn1-xn)'*(xn1-xn))) * (xn1-xn); 
                DG_T_x = 0.5*(DG_1_T_x_pn + DG_1_T_x_pn1);

            DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(xn1 -xn)) / ((xn1-xn)'*(xn1-xn))) * (xn1-xn);

            else
                % use MP evaluation if qn1 is approx. qn
                DG_T_x = D_1_T_n05;
                DG_Vext = DVext_n05;
            end
            
            % discrete gradient of kinetic energy w.r.t velocity 
            DG_T_p = 0.5*(IMn + IMn1)*p_n05;
            DG_T_x = [DG_T_x;0;0];
            DG_Vext = [DG_Vext;0;0];


            %% Residual vector
            resi = [qn1 - qn - h * DG_T_p; 
                    pn1 - pn + h * DG_T_x + h * DG_Vext];

            %% Tangent matrix
            tang = [];

        end

    end

end