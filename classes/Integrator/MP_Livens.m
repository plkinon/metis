classdef MP_Livens < Integrator

    %% Midpoint-Integration scheme for standard constrained DAE
    %
    % - not derived from variational principle
    %
    % - uses Livens equations of motion
    %
    % - takes account of non-constant mass-matrices
    %
    % Author: Philipp Kinon
    % Date  : 12.01.2024

    methods

        function self = MP_Livens(this_simulation, this_system)
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

            %% Abbreviations
            h = self.DT;
            nDOF = this_system.nDOF;
            mConstraints = this_system.mCONSTRAINTS;
            %% Unknows which will be iterated
            qn1 = zn1(1:nDOF);
            pn1 = zn1(nDOF+1:2*nDOF);
            vn1 = zn1(2*nDOF+1:3*nDOF);
            lambdan1 = zn1(3*nDOF+1:3*nDOF+mConstraints);
            g_n1 = this_system.constraint(qn1);

            %% Known quantities from last time-step
            qn = zn(1:nDOF);
            pn = zn(nDOF+1:2*nDOF);
            vn = zn(2*nDOF+1:3*nDOF);

            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            p_n05 = 0.5 * (pn + pn1);
            v_n05 = 0.5 * (vn + vn1);
            DVext_n05 = this_system.external_potential_gradient(q_n05);
            Mn05 = this_system.get_mass_matrix(q_n05);
            DVint_n05 = this_system.internal_potential_gradient(q_n05);
            G_n05 = this_system.constraint_gradient(q_n05);
            DT_q_n05 = this_system.kinetic_energy_gradient_from_velocity(q_n05, v_n05);

            %% Residual vector
             resi = [qn1 - qn - h * v_n05; 
                     pn1 - pn + h * DVext_n05 + h * DVint_n05 - h * DT_q_n05 + h * G_n05' * lambdan1; 
                     p_n05 - Mn05*v_n05;
                     g_n1];
            %% Tangent matrix
            tang = [];

        end

    end

end