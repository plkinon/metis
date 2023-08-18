classdef EM_pH < Integrator
    % Energy_Momentum-Integration scheme for standard constrained DAE
    %
    % - not derived from variational principle
    %
    % - uses discrete gradient for ext. potential and internal potential
    %
    % - uses mixed variables (e.g. for strains) and has port-Hamiltonian
    %   structure
    %
    % Author: Philipp Kinon
    % Date  : 26.04.2023

    methods

        function self = EM_pH(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF + this_system.mMixedQuantities;
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'EM-PH';
            self.has_enhanced_constraint_force = [];
            self.compute_potential_from_mixed_quantity = true;
end

        function z0 = set_initial_condition(~, this_simulation, this_system)
            % sets initial conditions

            z0 = [this_simulation.Q_0', (this_system.get_mass_matrix(this_simulation.Q_0) * this_simulation.V_0)', this_simulation.V_0', this_simulation.ALPHA_0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)
            % computes residual tangent

            %% Abbreviations
            h = self.DT;
            n = this_system.nDOF;
            p = this_system.nPotentialInvariants;
            mMixed = this_system.mMixedQuantities;
            M = this_system.MASS_MAT;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);
            Cn1 = zn1(3*n+1:3*n+mMixed);
            Vext_n1 = this_system.external_potential(qn1);
            Vint_n1 = this_system.internal_potential_from_mixed_quantity(Cn1);
            
            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            vn = zn(2*n+1:3*n);
            Cn = zn(3*n+1:3*n+mMixed);
            Vext_n = this_system.external_potential(qn);
            Vint_n = this_system.internal_potential_from_mixed_quantity(Cn);
           
            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            v_n05 = 0.5 * (vn + vn1);
            C_n05 = 0.5 * (Cn + Cn1);

            DVext_n05 = this_system.external_potential_gradient(q_n05);
            DVint_n05 = this_system.internal_potential_gradient_from_mixed_quantity(C_n05);

            %% Discrete gradients
            
            % for the internal potential
            if abs((Cn1-Cn)'*(Cn1-Cn)) > 1e-10
                DG_Vint = (Vint_n1 - Vint_n)/(Cn1-Cn);
            else
                % else use MP evaluation of gradient
                DG_Vint = DVint_n05;
            end

            % for the external potential
            if abs((qn1-qn)'*(qn1-qn)) > 1e-10               
                DG_Vext = DVext_n05 + ((Vext_n1 - Vext_n - DVext_n05'*(qn1 -qn)) / ((qn1-qn)'*(qn1-qn)) ) * (qn1-qn);
            else
                DG_Vext = DVext_n05;
            end
            D_C_q_n05 = this_system.mixed_quantity_gradient(q_n05);

%             if any(this_system.isCyclicCoordinate) 
%                 DG_1_T(this_system.isCyclicCoordinate) = zeros(size(DG_1_T(this_system.isCyclicCoordinate)));
%                 DG_Vext(this_system.isCyclicCoordinate) = zeros(size(DG_Vext(this_system.isCyclicCoordinate)));
%                 DG_Vint(this_system.isCyclicCoordinate) = zeros(size(DG_Vint(this_system.isCyclicCoordinate)));
%             end

            %% Residual vector
            resi = [qn1 - qn - h * v_n05; 
                    M*vn1 - M*vn + h * DG_Vext + h * D_C_q_n05 * DG_Vint;
                    Cn1 - Cn - h *  D_C_q_n05' * v_n05;
                    pn1 - M*vn1];

            %% Tangent matrix
            tang = [];
        end

    end

end