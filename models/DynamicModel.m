classdef DynamicModel < handle
    % DYNAMICMODEL Abstract interface for continuous-time physics models.
    
    properties
        deltaT      % Sampling time
        n           % Number of states
        m           % Number of inputs
    end
    
    methods (Abstract)
        % Returns state derivative x_dot given state x and input u
        dxdt = continuousDynamics(obj, x, u)
        
        % Computes analytical discrete-time linearization (Ad, Bd, cd) around xbar, ubar
        [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
    end
    
    methods
        function x_next = discreteStep(obj, x, u)
            % discreteStep: Standard RK4 forward integration
            dt = obj.deltaT;
            
            k1 = obj.continuousDynamics(x, u);
            k2 = obj.continuousDynamics(x + k1 * (dt / 2), u);
            k3 = obj.continuousDynamics(x + k2 * (dt / 2), u);
            k4 = obj.continuousDynamics(x + k3 * dt, u);
            
            x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
        end
    end
end
