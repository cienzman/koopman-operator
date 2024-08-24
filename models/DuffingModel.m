classdef DuffingModel < DynamicModel
    % DUFFINGMODEL Duffing Oscillator exhibiting double-well potential and
    % potentially chaotic dynamics.
    
    properties
        delta       % Damping parameter
        alpha       % Linear spring stiffness (negative for double-well)
        beta        % Nonlinear spring stiffness
    end
    
    methods
        function obj = DuffingModel(deltaT, delta_p, alpha_p, beta_p)
            obj.deltaT = deltaT;
            obj.n = 2;
            obj.m = 1;
            
            if nargin > 1
                obj.delta = delta_p;
                obj.alpha = alpha_p;
                obj.beta = beta_p;
            else
                % Default to standard double-well settings
                obj.delta = 0.2;
                obj.alpha = -1;
                obj.beta = 1;
            end
        end
        
        function dxdt = continuousDynamics(obj, x, u)
            % continuousDynamics: evaluates x_dot = f(x, u)
            % x: [2 x N] state matrix
            dxdt = zeros(2, size(x, 2));
            dxdt(1, :) = x(2, :);
            dxdt(2, :) = -obj.delta * x(2, :) - obj.alpha * x(1, :) - obj.beta * (x(1, :).^3) + u;
        end
        
        function [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
            % Continuous Jacobians evaluated at (xbar, ubar)
            Ac = [0, 1;
                 -obj.alpha - 3*obj.beta*(xbar(1)^2), -obj.delta];
            Bc = [0; 1];
            
            c_val = obj.continuousDynamics(xbar, ubar);
            
            % Zero-Order Hold exact discretization via matrix exponential
            M = [Ac, Bc, c_val; 
                 zeros(2, 4)];
            ABc_d = expm(M * obj.deltaT);
            
            Ad = ABc_d(1:2, 1:2);
            Bd = ABc_d(1:2, 3);
            cd = ABc_d(1:2, 4);
        end
    end
end
