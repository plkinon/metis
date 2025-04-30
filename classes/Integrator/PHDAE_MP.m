classdef PHDAE_MP < Integrator
    % Midpoint Integration scheme for port-Hamiltonian
    % semiexplicit DAEs
    %
    % - based only on constraint on velocity level
    %
    % - not derived from variational principle
    %
    % - assumes constant mass matrix
    
    methods
        
        function self = PHDAE_MP(this_simulation, this_system)
            self.DT = this_simulation.DT;
            self.T_0 = this_simulation.T_0;
            self.T_END = this_simulation.T_END;
            self.t = this_simulation.T_0:this_simulation.DT:this_simulation.T_END;
            self.NT = size(self.t, 2) - 1;
            self.nVARS = 2 * this_system.nDOF + 1 * this_system.mCONSTRAINTS;
            self.INDI_VELO = true;
            self.LM0 = zeros(this_system.mCONSTRAINTS, 1);
            self.hasPARA = false;
            self.NAME = 'PHDAE-MP';
            self.has_enhanced_constraint_force = false;
        end
        
        function z0 = set_initial_condition(self, this_simulation, ~)
            q0 = this_simulation.Q_0;
            v0 = this_simulation.V_0;
            z0 = [q0', v0', self.LM0'];
            
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
            nDOF = this_system.nDOF;
            m = this_system.mCONSTRAINTS;
            p = this_system.nPotentialInvariants;
            
            %% Unknows which will be iterated
            qn1 = zn1(1:nDOF);
            vn1 = zn1(nDOF+1:2*nDOF);
            lambda_n1 = zn1(2*nDOF+1:2*nDOF+m);
            M = this_system.get_mass_matrix(qn1);
            
            %% Known quantities from last time-step
            qn = zn(1:nDOF);
            vn = zn(nDOF+1:2*nDOF);
            
            %% MP evaluated quantities
            q_n05 = 0.5 * (qn + qn1);
            v_n05 = 0.5 * (vn + vn1);
            DVext_n05 = this_system.external_potential_gradient(q_n05);
            DVint_n05 = this_system.internal_potential_gradient(q_n05);
            G_n05 = this_system.constraint_gradient(q_n05);
            
            if ismethod(this_system,'get_dissipation_matrix')
                Dn05 = this_system.get_dissipation_matrix(q_n05);
            end

            %% Residual vector
            resi = [qn1 - qn - h * v_n05;
                M*(vn1 - vn) + h * DVext_n05 + h * DVint_n05 + h * Dn05 * v_n05 + h * G_n05' * lambda_n1;
                G_n05*v_n05];
            
            %% Tangent matrix
            tang = [];
            
        end
        
    end
    
end