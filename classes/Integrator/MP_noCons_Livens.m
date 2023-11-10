classdef MP_noCons_Livens < Integrator
    % Midpoint-Integration scheme for standard constrained DAE
    %
    % - based only on constraint on position level
    %
    % - independent momenta variables (Livens approach)
    %
    % - not derived from variational principle but simply evaluates RHS at
    %   t_{n+1/2}
    %


    methods

        function self = MP_noCons_Livens(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF;
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'MP-noCons';
            self.has_enhanced_constraint_force = [];
        end

        function z0 = set_initial_condition(~, this_simulation, this_system)

            z0 = [this_simulation.Q_0', (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0)',  this_simulation.V_0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)
            % Computes residual vector & tangent matrix
            %
            % :param zn1: state vector for next time step
            % :param zn: state vector at current time step
            % :param this_system: System object
            % :returns: [ResidualVector, TangentMatrix] for the Newton's method to update zn1

            %% Abbreviations
            
            h = self.DT;
            n = this_system.nDOF;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            vn = zn(2*n+1:3*n);

            %% Alpha evaluated quantities
            q_n05 = 0.5 * qn + 0.5 * qn1;
            p_n05 = 0.5 * pn + 0.5 * pn1;
            v_n05 = 0.5 * vn + 0.5 * vn1;
            Mn05 = this_system.get_mass_matrix(q_n05);
            DV_n05 = this_system.internal_potential_gradient(q_n05) + this_system.external_potential_gradient(q_n05);
            D2V_n05 = this_system.internal_potential_hessian(q_n05) + this_system.external_potential_hessian(q_n05);
            Dq_T_n05 = this_system.kinetic_energy_gradient_from_velocity(q_n05, v_n05);
            %% Residual vector
            resi = [qn1 - qn - h * v_n05; 
                    pn1 - pn - h * Dq_T_n05 + h * DV_n05;
                    p_n05 - Mn05*v_n05];

            %% Tangent matrix
            %tang = [eye(n), -h * 0.5 * IMn05; 
            %       h * 0.5 * D2V_n05 , eye(n)];
            tang = [];
        end

    end

end