classdef EML_noCons_Lagrange < Integrator
    % Energy_Momentum-Integration scheme for standard constrained DAE
    %
    % - not derived from variational principle
    %
    % - uses Livens equations of motion
    %
    % - uses discrete gradient for Lagrangians
    %
    % - takes account of non-constant mass-matrices
    %
    % Author: Philipp Kinon
    % Date  : 09.12.2022

    methods

        function self = EML_noCons_Lagrange(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF;
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'EML-noCons-Lagrange';
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

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);
            Mn1 = this_system.get_mass_matrix(qn1);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            vn = zn(2*n+1:3*n);
            Mn = this_system.get_mass_matrix(qn);

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            v_n05 = 0.5 * (vn + vn1);
            Mn05 = this_system.get_mass_matrix(q_n05);

            D_1_L_n05 = this_system.kinetic_energy_gradient_from_velocity(q_n05, v_n05) - this_system.external_potential_gradient(q_n05) - this_system.internal_potential_gradient(q_n05);
            
            L_qn1vn = 1/2 * vn' * this_system.get_mass_matrix(qn1) * vn - this_system.external_potential(qn1) - this_system.internal_potential(qn1); %this_system.lagrange_function(qn1, vn);
            L_qnvn  = 1/2 * vn' * this_system.get_mass_matrix(qn) * vn - this_system.external_potential(qn) - this_system.internal_potential(qn); %this_system.lagrange_function(qn, vn);
            L_qn1vn1 = 1/2 * vn1' * this_system.get_mass_matrix(qn1) * vn1 - this_system.external_potential(qn1) - this_system.internal_potential(qn1); %this_system.lagrange_function(qn1, vn1);
            L_qnvn1 = 1/2 * vn1' * this_system.get_mass_matrix(qn) * vn1 - this_system.external_potential(qn) - this_system.internal_potential(qn); %this_system.lagrange_function(qn, vn1);
            
            D_1_L_qn05_vn = this_system.kinetic_energy_gradient_from_velocity(q_n05, vn) - this_system.external_potential_gradient(q_n05) - this_system.internal_potential_gradient(q_n05);
            D_1_L_qn05_vn1 = this_system.kinetic_energy_gradient_from_velocity(q_n05, vn1) - this_system.external_potential_gradient(q_n05) - this_system.internal_potential_gradient(q_n05);

            if abs((qn1-qn)'*(qn1-qn)) > 1e-9
                DG_1_L_q_vn = D_1_L_qn05_vn + ((L_qn1vn - L_qnvn - D_1_L_qn05_vn'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_L_q_vn1 = D_1_L_qn05_vn1 + ((L_qn1vn1 - L_qnvn1 - D_1_L_qn05_vn1'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn))) * (qn1-qn); 
                DG_1_L = 0.5*(DG_1_L_q_vn + DG_1_L_q_vn1);
                DG_2_L = 0.5*(Mn + Mn1)*v_n05;
            else
                DG_1_L = D_1_L_n05;
                DG_2_L = Mn05*v_n05;
            end

%             if any(this_system.isCyclicCoordinate) 
%                 DG_1_L(this_system.isCyclicCoordinate) = D_1_L_n05(this_system.isCyclicCoordinate);
%             end

            %% Residual vector
            resi = [qn1 - qn - h * v_n05; 
                    pn1 - pn - h * DG_1_L;
                    p_n05 - DG_2_L];

            %% Tangent matrix
            tang = [];
        end

    end

end