classdef LinearImplicit_pH < Integrator
    % Energy_Momentum-Integration scheme for PH mechanical system
    %
    % - not derived from variational principle
    %
    % - uses discrete gradient for ext. potential and internal potential
    %
    % - uses mixed variables (e.g. for strains) and has port-Hamiltonian
    %   structure
    %
    % - more info: https://doi.org/10.1002/pamm.202300144

    methods

        function self = LinearImplicit_pH(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF + this_system.mMixedQuantities;
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'LinearImplicit-PH';
            self.has_enhanced_constraint_force = [];
            self.compute_potential_from_mixed_quantity = true;
end

        function z0 = set_initial_condition(self, this_simulation, this_system)

            z0 = [this_simulation.Q_0'-this_simulation.V_0'*self.DT/2, (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0)', this_simulation.V_0', this_simulation.ALPHA_0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system, ~)
            % Computes residual vector & tangent matrix
            %
            % :param zn1: state vector for next time step
            % :param zn: state vector at current time step
            % :param this_system: System object
            % :returns: [ResidualVector, TangentMatrix] for the Newton's method to update zn1

            %% Abbreviations
            h = self.DT;
            n = this_system.nDOF;
            mMixed = this_system.mMixedQuantities;
            M = this_system.MASS_MAT;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);
            Cn1 = zn1(3*n+1:3*n+mMixed);
            
            %% Known quantities from last time-step
            qn = zn(1:n);
            vn = zn(2*n+1:3*n);
            Cn = zn(3*n+1:3*n+mMixed);
           
            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            v_n05 = 0.5 * (vn + vn1);
            C_n05 = 0.5 * (Cn + Cn1);
            q_n_bar = qn + h*vn;

            DVext_n05 = this_system.external_potential_gradient(q_n05);
            DVint_n05 = this_system.internal_potential_gradient_from_mixed_quantity(C_n05);
            D_C_q_n_bar = this_system.mixed_quantity_gradient(q_n_bar);
            D2Vint_n05 = this_system.internal_potential_hessian_from_mixed_quantity(C_n05);
            D_diss_bar = this_system.get_dissipation_matrix(q_n_bar);

            %% Residual vector
            resi = [qn1 - qn - h * vn; 
                    M*vn1 - M*vn + h * DVext_n05 + h * D_C_q_n_bar * DVint_n05 + h*D_diss_bar*v_n05;
                    Cn1 - Cn - h *  D_C_q_n_bar' * v_n05;
                    pn1 - M*vn1];
            
            lenq = size(qn1,1);
            lenC = size(Cn1,1);

            %% Tangent matrix
            tang = [eye(lenq,lenq), zeros(lenq,lenq), zeros(lenq,lenq), zeros(lenq,lenC); %check
                    zeros(lenq,lenq), zeros(lenq,lenq), M+h/2*D_diss_bar, h * D_C_q_n_bar*D2Vint_n05*1/2; %check
                    zeros(lenC,lenq), zeros(lenC,lenq), -h*D_C_q_n_bar'*1/2, eye(lenC,lenC);
                    zeros(lenq,lenq), eye(lenq,lenq),   -M, zeros(lenq, lenC)          ];
        end

    end

end
