classdef Integrator
    % Abstract Integrator class

    properties
        NAME
        DT
        T_0
        T_END
        t
        NT
        nVARS
        INDI_VELO
        LM0
        hasPARA
        PARA
        has_enhanced_constraint_force
        compute_potential_from_mixed_quantity
    end

end
