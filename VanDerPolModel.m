classdef VanDerPolModel
    % VANDERPOLMODEL Encapsulates the forced Van der Pol oscillator.
    %
    %   Continuous dynamics:
    %       x1_dot = 2 * x2
    %       x2_dot = -0.8*x1 + 2*x2 - 10*x1^2*x2 + u
    %
    %   Discrete-time stepping uses 4th-order Runge-Kutta (RK4).
    %   Local linearizations use zero-order-hold (ZOH) exact discretization.

    properties
        deltaT          % Sampling period [s]
        n = 2;          % State dimension
        m = 1;          % Input dimension
    end

    methods
        function obj = VanDerPolModel(deltaT)
            obj.deltaT = deltaT;
        end

        function dxdt = continuousDynamics(~, x, u)
            % Evaluates the vector field  f(x, u)  at N points simultaneously.
            % x : [2 x N],  u : [1 x N]
            dxdt = [2 .* x(2, :);
                   -0.8 .* x(1, :) + 2 .* x(2, :) - 10 .* x(1,:).^2 .* x(2,:) + u];
        end

        function x_next = discreteStep(obj, x, u)
            % RK4 integration over one sampling period.
            dt = obj.deltaT;
            k1 = obj.continuousDynamics(x,              u);
            k2 = obj.continuousDynamics(x + dt/2 * k1,  u);
            k3 = obj.continuousDynamics(x + dt/2 * k2,  u);
            k4 = obj.continuousDynamics(x + dt   * k3,  u);
            x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        end

        function [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
            % Returns the affine discrete-time model linearized at (xbar, ubar):
            %
            %   x_{k+1} = Ad*(x_k - xbar) + Bd*(u_k - ubar) + cd + xbar
            %
            % Computed via ZOH exact discretization of the continuous Jacobians.

            % Continuous Jacobian  df/dx  evaluated at (xbar, ubar)
            Ac = [0,                          2;
                 -0.8 - 20*xbar(1)*xbar(2),  2 - 10*xbar(1)^2];
            Bc = [0; 1];

            % Equilibrium drift at the linearization point
            fc = obj.continuousDynamics(xbar, ubar);

            % ZOH via matrix exponential of the augmented system [A B c; 0]
            Maug = [Ac, Bc, fc;
                    zeros(2, 4)];
            Md = expm(Maug * obj.deltaT);

            Ad = Md(1:2, 1:2);
            Bd = Md(1:2, 3);
            cd = Md(1:2, 4);
        end
    end
end