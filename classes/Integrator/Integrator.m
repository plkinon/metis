classdef Integrator
    % Abstract Integrator class

    properties
        NAME % integrator name
        DT % time step size 
        T_0 % start time
        T_END % end time
        t % time
        NT % number of time steps
        nVARS % number of unknown variables
        INDI_VELO % whether integrator deals with independent velocities AND momenta
        LM0 % initial value for Lagrange multipliers
        hasPARA % whether integrator features parameters
        PARA % integrator parameters
        has_enhanced_constraint_force % whether constraint force is enhances
        compute_potential_from_mixed_quantity % whether potential is computed from mixed quantity or not
    end

end
