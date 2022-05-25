classdef CSE_B < Integrator

    %% Variational integration scheme for constrained DAE
    %
    % - based on constraint on position level
    %
    % - generalization of the symplectic Euler B method
    %
    % - independent momentum variables
    %
    % - constraints are enforced at t_{n+1}
    %
    % Author: Philipp Kinon
    % Date  : 09.12.2020

    methods

        function self = CSE_B(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 2 * this_system.nDOF + 1 * this_system.mCONSTRAINTS;
            self.INDI_VELO = false;
            self.LM0 = zeros(1*this_system.mCONSTRAINTS, 1);
            self.hasPARA = false;
            self.NAME = 'CSE-B';
            self.has_enhanced_constraint_force = false;
        end

        function z0 = set_initial_condition(self, this_simulation, this_system)

            z0 = [this_simulation.Q_0', (this_system.MASS_MAT * this_simulation.V_0)', self.LM0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)

            %% Abbreviations
            M = this_system.MASS_MAT;
            IM = M \ eye(size(M));
            h = self.DT;
            n = this_system.nDOF;
            m = this_system.mCONSTRAINTS;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            lambdan = zn1(2*n+1:2*n+m);
            G_n1 = this_system.constraint_gradient(qn1);
            g_n1 = this_system.constraint(qn1);


            % Hessian of constraints are multiplied by inverse MassMat and
            % pn1 for each constraint to avoid 3rd order tensors
            T_n1 = zeros(m, n);
            for l = 1:m
                tmp = this_system.constraint_hessian(qn1, l);
                for k = 1:n
                    T_n1(l, k) = tmp(:, k)' * IM * pn1;
                end
            end

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            G_n = this_system.constraint_gradient(qn);
            DV_n = this_system.internal_potential_gradient(qn) + this_system.external_potential_gradient(qn);

            %% Residual vector
            resi = [qn1 - qn - h * IM * pn1; ...
                pn1 - pn + h * DV_n + h * G_n' * lambdan; ...
                g_n1];

            %% Tangent matrix
            tang = [eye(n), -h * IM, zeros(n, m); ...
                zeros(n), eye(n), h * G_n'; ...
                G_n1, zeros(n, m)', zeros(m)];

        end

    end

end