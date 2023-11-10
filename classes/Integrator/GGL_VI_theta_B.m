classdef GGL_VI_theta_B < Integrator
    % Runge-Kutta typed scheme for GGL-like constrained DAE
    %
    % - based on constraint on position and velocity level
    %   (GGL-stabilisation)
    %
    % - independent momentum and velocity variables (Hamilton Potryagin approach)
    %
    % - derived from variational principle
    %
    % - from symplectic-theta-framework 
    %
    % - more info: https://doi.org/10.1007/s11071-023-08522-7

    methods

        function self = GGL_VI_theta_B(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF + 2 * this_system.mCONSTRAINTS;
            self.INDI_VELO = true;
            self.LM0 = zeros(2*this_system.mCONSTRAINTS, 1);
            self.hasPARA = true;
            self.PARA = this_simulation.INT_PARA(:); %[The:round theta  theta: vartheta]; [0.5 0.5]: more stable, [1 0.5]: exact constraint vel. level
            self.NAME = 'GGL-VI-theta-B';
            self.has_enhanced_constraint_force = true;
        end

        function z0 = set_initial_condition(self, this_simulation, this_system)

            M = this_system.MASS_MAT;
            v0 = this_simulation.V_0;
            q0 = this_simulation.Q_0;           
            p0 = M * v0;
            z0 = [q0', p0', v0', self.LM0'];

        end

        function [resi, tang] = compute_resi_tang(self, zn1, zn, this_system)
            % Computes residual vector & tangent matrix
            %
            % :param zn1: state vector and next time step
            % :param zn: state vector at current time step
            % :param this_system: System object
            % :returns: [ResidualVector, TangentMatrix] for the Newton's method to update zn1

            %% Abbreviations
            M = this_system.MASS_MAT;
            IM = M \ eye(size(M));
            h = self.DT;
            n = this_system.nDOF;
            m = this_system.mCONSTRAINTS;

            %% Unknows which will be iterated
            qn1 = zn1(1:n);
            pn1 = zn1(n+1:2*n);
            vn1 = zn1(2*n+1:3*n);
            lambdan = zn1(3*n+1:3*n+m);
            gamman = zn1(3*n+m+1:end);
            G_n1 = this_system.constraint_gradient(qn1);
            g_n1 = this_system.constraint(qn1);

            % Hessian of constraints are multiplied by inverse MassMat and
            % pn1 for each constraint to avoid 3rd order tensors
            T_n1 = zeros(m, n);
            t_n1 = zeros(n);
            for l = 1:m
                tmp = this_system.constraint_hessian(qn1, l);
                t_n1 = t_n1 + this_system.constraint_hessian(qn1, l) * lambdan(l);
                for k = 1:n
                    T_n1(l, k) = tmp(:, k)' * IM * pn1;
                end
            end

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            G_n = this_system.constraint_gradient(qn);

            %% Quantities at t_n+theta
            theta = self.PARA(2);
            The = self.PARA(1);
            q_nt = (1 - The) * qn + The * qn1;
            p_n1mt = The * pn + (1 - The) * pn1;
            DV_nt = this_system.internal_potential_gradient(q_nt) + this_system.external_potential_gradient(q_nt);
            G_nt = this_system.constraint_gradient(q_nt);
            D2V_nt = this_system.internal_potential_hessian(q_nt) + this_system.external_potential_hessian(q_nt);

            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_nt = zeros(n);
            T_nt = zeros(m, n);
            for j = 1:m
                tmp2 = this_system.constraint_hessian(q_nt, j);
                t_nt = t_nt + this_system.constraint_hessian(q_nt, j) * gamman(j);
                for k = 1:n
                    T_nt(j, k) = tmp2(:, k)' * vn1;
                end
            end

            %% Residual vector
            resi = [qn1 - qn - h * vn1 - h * IM * G_nt' * gamman; pn1 - pn + h * DV_nt + h * ((1 - theta) * G_n' + theta * G_n1') * lambdan + h * t_nt * vn1; (M * vn1 - p_n1mt + h * (The * (1 - theta) * G_n' - theta * (1 - The) * G_n1') * lambdan); g_n1; G_nt * vn1];

            %% Tangent matrix
            tang = [eye(n) - h * IM * The * t_nt, zeros(n), -h * eye(n), zeros(n, m), -h * IM * G_nt'; h * D2V_nt * The + h * theta * t_n1, eye(n), h * t_nt, h * ((1 - theta) * G_n' + theta * G_n1'), h * T_nt'; h * theta * (1 - The) * t_n1, -(1 - The) * eye(n), M, h * (The * (1 - theta) * G_n' - theta * (1 - The) * G_n1'), zeros(n, m); G_n1, zeros(n, m)', zeros(n, m)', zeros(m), zeros(m); The * T_nt, zeros(n, m)', G_nt, zeros(m), zeros(m)];
            %tang = [];
        end

    end

end