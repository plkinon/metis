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
        LM0
    end


    methods

        function self = Integrator(CONFIG)
            self.DT    = CONFIG.DT;
            self.T_0   = CONFIG.T_0;
            self.T_END = CONFIG.T_END;
            self.t     = CONFIG.T_0:CONFIG.DT:CONFIG.T_END;
            self.NT    = size(self.t, 2) - 1;
        end
    end
end
