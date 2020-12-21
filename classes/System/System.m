classdef System
    % Abstract system class

    properties
        MASS
        MASS_MAT
        EXT_ACC
        DIM
        nDOF
        mCONSTRAINTS
        nBODIES
        nLM
        GEOM
        nPotentialInvariants
        nConstraintInvariants
        nVconstraintInvariants
        z
        H
        T
        V
        J
        Hdiff
        Jdiff
        constraint_position
        constraint_velocity
    end

end
