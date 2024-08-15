classdef VanDerPolModel
    % VANDERPOLMODEL Encapsulates the Van der Pol oscillator dynamics.
    % Provides methods for continuous-time evaluation, discrete-time 
    % Runge-Kutta 4 stepping, and local linearizations about a given state.
    
    properties
        deltaT      % Sampling time
        n = 2;      % Number of states
        m = 1;      % Number of inputs
    end
    
    methods
        function obj = VanDerPolModel(deltaT)
            if nargin < 1
                obj.deltaT = 0.01;
            else
                obj.deltaT = deltaT;
            end
        end
        
        function dxdt = continuousDynamics(~, x, u)
            % continuousDynamics: evaluates x_dot = f(x, u)
            % x: [2 x N] state matrix
            % u: [1 x N] input matrix
            dxdt = zeros(2, size(x, 2));
            dxdt(1, :) = 2 * x(2, :);
            dxdt(2, :) = -0.8 * x(1, :) + 2 * x(2, :) - 10 * (x(1, :).^2) .* x(2, :) + u;
        end
        
        function x_next = discreteStep(obj, x, u)
            % discreteStep: RK4 forward integration
            dt = obj.deltaT;
            
            k1 = obj.continuousDynamics(x, u);
            k2 = obj.continuousDynamics(x + k1 * (dt / 2), u);
            k3 = obj.continuousDynamics(x + k2 * (dt / 2), u);
            k4 = obj.continuousDynamics(x + k3 * dt, u);
            
            x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
        end
        
        function [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
            % localLinearization: Computes analytical discrete-time linearization
            % about an operating point (xbar, ubar).
            % Returns matrices for model: x_k+1 = Ad*(x_k - xbar) + Bd*(u_k - ubar) + cd + xbar
            
            % Continuous Jacobians evaluated at (xbar, ubar)
            Ac = [0, 2;
                 -0.8 - 20*xbar(1)*xbar(2), 2 - 10*xbar(1)^2];
            Bc = [0; 1];
            
            c_val = obj.continuousDynamics(xbar, ubar);
            
            % Zero-Order Hold exact discretization via matrix exponential
            % of the augmented matrix
            M = [Ac, Bc, c_val; 
                 zeros(2, 4)];
            ABc_d = expm(M * obj.deltaT);
            
            Ad = ABc_d(1:2, 1:2);
            Bd = ABc_d(1:2, 3);
            cd = ABc_d(1:2, 4);
        end
    end
end
