classdef Lagrange_top_ODE < Integrator
    % Integration scheme for Lagrange top 
    %
    % - based on angular momentum and director d3 (ODE approach)
    %
    % - not derived from variational principle 
    %
    % - used for NODY publication, taken from Bobenko & Suris, 1999
    %
    % Author: Philipp Kinon
    % Date  : 24.01.2023

    methods

        function self = Lagrange_top_ODE(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 2*this_system.nDOF; % 3 components of d_3 and 3 components of m (angular momentum w.r.t. origin)
            self.INDI_VELO = true;
            self.LM0 = [];
            self.hasPARA = false;
            self.NAME = 'Lagrange_top_ODE';
            self.has_enhanced_constraint_force = [];
        end

        function z0 = set_initial_condition(~, this_simulation, ~)

            z0 = [this_simulation.Q_0', this_simulation.V_0'];

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
            g = abs(this_system.EXT_ACC(3));
            m = this_system.MASS;
            %% Unknows which will be iterated
            m_n1 = zn1(4:6);
            d3_n1 = zn1(1:3);

            %% Known quantities from last time-step
            m_n = zn(4:6);
            d3_n = zn(1:3);

            %% Alpha evaluated quantities
            d3_n05 = 0.5 * d3_n + 0.5 * d3_n1;
            
            %% Geometric parameters
            a = 0.1; % length of the gyro top
            r = a / 2; % radius of the gyro top
            l = 3 * a / 4; % location of center of mass along sym. axis
            J1 = 3 / 80 * m * (4 * r^2 + a^2); % inertia moment w.r.t. d1-axis (J1 = J2)
            % momenta w.r.t. origin in d_i-system
            I1 = J1 + m*l^2;
            I1_inv = 1/I1;

            %% Residual vector
            resi = [m_n1 - m_n - h * m*g*l * cross([0;0;1], d3_n); 
                    d3_n1 - d3_n - h * I1_inv * cross(m_n1, d3_n05)];

            %% Tangent matrix
            tang = [];
            
        end

    end

end