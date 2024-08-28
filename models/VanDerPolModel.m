classdef VanDerPolModel < DynamicModel
    % VANDERPOLMODEL Encapsulates the Van der Pol oscillator dynamics.
    
    methods
        function obj = VanDerPolModel(deltaT)
            if nargin < 1
                obj.deltaT = 0.01;
            else
                obj.deltaT = deltaT;
            end
            obj.n = 2;
            obj.m = 1;
        end
        
        function dxdt = continuousDynamics(obj, x, u)
            % continuousDynamics: evaluates x_dot = f(x, u)
            dxdt = zeros(2, size(x, 2));
            dxdt(1, :) = 2 * x(2, :);
            dxdt(2, :) = -0.8 * x(1, :) + 2 * x(2, :) - 10 * (x(1, :).^2) .* x(2, :) + u;
        end
        
        function [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
            % Continuous Jacobians evaluated at (xbar, ubar)
            Ac = [0, 2;
                 -0.8 - 20*xbar(1)*xbar(2), 2 - 10*xbar(1)^2];
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
