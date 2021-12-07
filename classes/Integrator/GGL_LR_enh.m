classdef GGL_LR_enh < Integrator

    %% Variational integration scheme for GGL-like constrained DAE
    %
    % - based on constraint on position and velocity level
    %   (GGL-stabilisation)
    %
    % - independent momenta variables (Hamilton Potryagin approach)
    %
    % - taken from Leimkuhler & Reich , p.189, using Gen-Alpha for Psi_h and
    %   making GGL-stabilisation (not a VI)
    %
    % - constraints are enforced at t_{n+1}
    %
    % Author: Philipp Kinon
    % Date  : 09.12.2020

    methods

        function self = GGL_LR_enh(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 4 * this_system.nDOF + 3 * this_system.mCONSTRAINTS;
            self.INDI_VELO = false;
            self.LM0 = zeros(3*this_system.mCONSTRAINTS, 1);
            self.hasPARA = true;
            self.PARA = this_simulation.INT_PARA(1);
            self.NAME = 'GGL-LR-enh';
        end

        function z0 = set_initial_condition(self, this_simulation, this_system)

            n = this_system.nDOF;
            z0 = [this_simulation.Q_0', (this_system.MASS_MAT * this_simulation.V_0)', zeros(n, 1)', zeros(n, 1)', self.LM0'];

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
            p_bar_n1 = zn1(2*n+1:3*n);
            p_bar_n = zn1(3*n+1:4*n);
            lambdan = zn1(4*n+1:4*n+m);
            gamma_n1 = zn1(4*n+m+1:4*n+2*m);
            kappa = zn1(4*n+2*m+1:end);
            G_n1 = this_system.constraint_gradient(qn1);
            g_n1 = this_system.constraint(qn1);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);

            %% Alpha-eval. quantities
            al = self.PARA;
            q_nal = al * qn1 + (1 - al) * qn;
            p_bar_n1mal = (1 - al) * p_bar_n1 + al * p_bar_n;
            G_n = this_system.constraint_gradient(qn);
            G_nal = this_system.constraint_gradient(q_nal);
            DV_nal = this_system.internal_potential_gradient(q_nal) + this_system.external_potential_gradient(q_nal);

            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n = zeros(n);
            for j = 1:m
                t_n = t_n + this_system.constraint_hessian(q_nal, j) * kappa(j);
            end

            %% Residual vector
            resi = [(qn1 - qn - h * IM * p_bar_n1mal - h * IM * G_nal' * kappa); (p_bar_n1 - p_bar_n + h * DV_nal + h * t_n * IM * p_bar_n1mal); (p_bar_n - pn + h / 2 * G_n' * lambdan); (p_bar_n1 - pn1 + h / 2 * G_n1' * gamma_n1); g_n1; G_n1 * IM * pn1; G_nal * IM * p_bar_n1mal];

            %% Tangent matrix
            tang = [];

        end

    end

end