classdef PaperVanDerPolModel < DynamicModel
    % PAPER VANDERPOL MODEL
    % Exactly replicates the Euler discretization from the paper:
    % "Stability of data-driven Koopman MPC with terminal conditions"
    
    properties
        nu = 0.1;
    end
    
    methods
        function obj = PaperVanDerPolModel()
            obj.deltaT = 0.05;
            obj.n = 2;
            obj.m = 1;
        end
        
        function dxdt = continuousDynamics(obj, x, u)
            % This model is explicitly defined in discrete time in the paper.
            % But we can provide the continuous drift for completeness.
            dxdt = zeros(2, size(x, 2));
            dxdt(1, :) = x(2, :);
            dxdt(2, :) = obj.nu * (1 - x(1, :).^2) .* x(2, :) - x(1, :) + u;
        end
        
        function x_next = discreteStep(obj, x, u)
            % Euler discretization as defined in the paper
            % x^+ = x + \Delta t * f_c(x, u)
            dxdt = obj.continuousDynamics(x, u);
            x_next = x + obj.deltaT * dxdt;
        end
        
        function [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
            % Numerical Jacobian for discrete step around (xbar, ubar)
            eps = 1e-6;
            Ad = zeros(obj.n, obj.n);
            Bd = zeros(obj.n, obj.m);
            
            for i = 1:obj.n
                x_eps = xbar;
                x_eps(i) = x_eps(i) + eps;
                Ad(:, i) = (obj.discreteStep(x_eps, ubar) - obj.discreteStep(xbar, ubar)) / eps;
            end
            
            for j = 1:obj.m
                u_eps = ubar;
                u_eps(j) = u_eps(j) + eps;
                Bd(:, j) = (obj.discreteStep(xbar, u_eps) - obj.discreteStep(xbar, ubar)) / eps;
            end
            
            cd = obj.discreteStep(xbar, ubar) - Ad * xbar - Bd * ubar;
        end
    end
end
