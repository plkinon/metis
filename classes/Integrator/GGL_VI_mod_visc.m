classdef GGL_VI_mod_visc < Integrator
    % Variational integration scheme for GGL-like constrained DAE
    %
    % - based on constraint on position and momentum level
    %
    % - independent momenta variables (Livens approach)
    %
    % - derived from variational principle 
    %
    % - symplectic
    %   
    % - constraints are enforced at t_{n+1}
    %
    % - takes into account non-conservative viscous forces
    %
    % - more info: https://doi.org/10.1007/s11044-023-09889-6
    methods

        function self = GGL_VI_mod_visc(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 3 * this_system.nDOF + 2 * this_system.mCONSTRAINTS;
            self.INDI_VELO = true;
            self.LM0 = zeros(2*this_system.mCONSTRAINTS, 1);
            self.hasPARA = false;
            self.NAME = 'GGL-VI (modified) ';
            self.has_enhanced_constraint_force = true;
        end

        function z0 = set_initial_condition(self, this_simulation, this_system)

            M = this_system.MASS_MAT;
            z0 = [this_simulation.Q_0', (M * this_simulation.V_0)', this_simulation.V_0', self.LM0'];

        end

        function z_rearranged = rearrange_unknowns(~, this_simulation, this_system)

            % v_n is an unknown of this scheme, has to be shifted backwards
            % by 1 after computation
            n = this_system.nDOF;
            z_rearranged = this_simulation.z;
            z_rearranged(1:(end -1), 2*n+1:3*n) = this_simulation.z(2:end, 2*n+1:3*n);
            z_rearranged(end, 2*n+1:3*n) = NaN;

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
            vn = zn1(2*n+1:3*n);
            lambdan = zn1(3*n+1:3*n+m);
            gamman1 = zn1(3*n+m+1:end);
            g_n1 = this_system.constraint(qn1);
            G_n1 = this_system.constraint_gradient(qn1);

            %% Known quantities from last time-step
            qn = zn(1:n);
            pn = zn(n+1:2*n);
            G_n = this_system.constraint_gradient(qn);
            DV_n = this_system.internal_potential_gradient(qn) + this_system.external_potential_gradient(qn);
            q_bar = qn + h * vn;
            G_bar = this_system.constraint_gradient(q_bar);
            F_nc_n = this_system.viscous_forces(vn);
            DF_nc_dv = this_system.viscous_forces_gradient(vn);

            % Hessian of constraints are multiplied by LMs for each
            % Constraint (avoid 3rd order tensor)
            t_n_bar = zeros(n);
            for j = 1:m
                t_n_bar = t_n_bar + this_system.constraint_hessian(q_bar, j) * gamman1(j);
            end
            T_bar = zeros(m, n);
            for l = 1:m
                tmp = this_system.constraint_hessian(q_bar, l);
                for k = 1:n
                    T_bar(l, k) = tmp(:, k)' * IM * pn1;
                end
            end

            %% Residual vector
            resi = [qn1 - qn - h * vn - h * IM * G_bar' * gamman1; 
                    pn1 - pn + h * DV_n + h * G_n' * lambdan + h * t_n_bar * IM * pn1 - h * F_nc_n; 
                    M * vn - pn1 - h * t_n_bar * IM * pn1; 
                    g_n1; 
                    G_bar * IM * pn1];

            %% Tangent matrix
            tang = [[eye(n), zeros(n), -h * eye(n) - h^2 * IM * t_n_bar, zeros(n, m), -h * IM * G_bar']; 
                    [zeros(n), eye(n) + h * t_n_bar * IM, -h*DF_nc_dv , h * G_n', h * T_bar']; 
                    [zeros(n), -eye(n) - h * t_n_bar * IM, M, zeros(n, m), -h * T_bar']; 
                    [G_n1, zeros(n, m)', zeros(n, m)', zeros(m), zeros(m)]; 
                    [zeros(n, m)', G_bar * IM, T_bar * h, zeros(m), zeros(m)]];


        end

    end

end