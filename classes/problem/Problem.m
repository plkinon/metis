classdef Problem
    % Abstract problem class

    properties
        MASS_MAT
        EXT_ACC
        DIM
        nDOF
        mCONSTRAINTS
        nLM
        GEOM
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
