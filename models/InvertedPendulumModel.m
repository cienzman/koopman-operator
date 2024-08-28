classdef InvertedPendulumModel < DynamicModel
    % INVERTEDPENDULUMMODEL Rotary Inverted Pendulum model.
    % State: [arm_angle, pendulum_angle, arm_velocity, pendulum_velocity]^T
    
    properties
        mb = 0.257;
        mp = 0.127;
        Lb = 0.216;
        Lp = 0.337;
        b1 = 0.009;
        b2 = 5.2e-4;
        g = 9.81;
        Rm = 2.6;
        kt = 0.00767;
        km = 0.00767;
        Kgi = 14;
        etag = 0.90;
        etam = 0.69;
    end
    
    methods
        function obj = InvertedPendulumModel(deltaT)
            if nargin < 1
                obj.deltaT = 0.01;
            else
                obj.deltaT = deltaT;
            end
            obj.n = 4;
            obj.m = 1;
        end
        
        function dxdt = continuousDynamics(obj, x, u)
            % continuousDynamics: Highly non-linear equations of motion
            % x: [4 x N] state matrix
            
            dxdt = zeros(4, size(x, 2));
            
            x2 = x(2, :);
            x3 = x(3, :);
            x4 = x(4, :);
            
            % Common terms
            cos_x2 = cos(x2);
            sin_x2 = sin(x2);
            cos_x2_sq = cos_x2.^2;
            
            alfa = 6*obj.Lp^3*obj.Rm*obj.mp + 8*obj.Lb^2*obj.Lp*obj.Rm*obj.mb + ...
                   24*obj.Lb^2*obj.Lp*obj.Rm*obj.mp - 6*obj.Lp^3*obj.Rm*obj.mp.*cos_x2_sq - ...
                   18*obj.Lb^2*obj.Lp*obj.Rm*obj.mp.*cos_x2_sq;
                   
            beta = 4*obj.Lb^2*obj.mb + 12*obj.Lb^2*obj.mp - ...
                   3*obj.Lp^2*obj.mp.*cos_x2_sq + 3*obj.Lp^2*obj.mp;
                   
            gamma = obj.Lp^2*obj.mp - ((9*obj.Lb^2*obj.Lp^2*obj.mp^2.*cos_x2_sq) ./ beta);
            
            % f3 calculations (arm acceleration)
            f3_1 = (9*obj.Lb*obj.Lp^2*obj.mp*obj.Rm.*(sin_x2 - sin_x2.^3)) .* (x3.^2);
            f3_2 = (6*obj.Lp^3*obj.mp*obj.Rm.*sin(2*x2)) .* x3 .* x4;
            f3_3 = (-24*obj.b1*obj.Lp*obj.Rm - 24*obj.Kgi^2*obj.km*obj.kt*obj.Lp*obj.etag*obj.etam) .* x3;
            f3_4 = (12*obj.Lb*obj.Lp^2*obj.mp*obj.Rm.*sin_x2) .* (x4.^2);
            f3_5 = (36*obj.b2*obj.Lb*obj.Rm.*cos_x2) .* x4;
            f3_6 = (24*obj.etag*obj.etam*obj.Kgi*obj.kt*obj.Lp.*u + 9*obj.g*obj.Lb*obj.Lp*obj.mp*obj.Rm.*sin(2*x2));
            
            f3 = (f3_1 + f3_2 + f3_3 + f3_4 + f3_5 + f3_6) ./ alfa;
            
            % f4 calculations (pendulum acceleration)
            f4_1 = 3*obj.Lp^2*obj.mp.*cos_x2.*sin_x2.*beta .* (x3.^2);
            f4_2 = -36*obj.Lb*obj.Lp^3*obj.mp^2.*cos_x2_sq.*sin_x2 .* x3 .* x4;
            f4_3 = -82*obj.Lb*obj.Lp*obj.mp.*cos_x2 .* (obj.b1 + (obj.etag*obj.etam*obj.Kgi^2*obj.km*obj.kt)/obj.Rm) .* x3;
            f4_4 = -36*obj.Lb^2*obj.Lp^2*obj.mp^2.*sin_x2.*cos_x2 .* (x4.^2);
            f4_5 = -12*obj.b2.*beta.*x4;
            f4_6 = ((6*obj.etag*obj.etam*obj.Kgi*obj.kt*obj.Lb*obj.Lp*obj.mp.*cos_x2.*u) ./ (obj.Rm.*beta)) + ...
                   0.5*obj.g*obj.Lp*obj.mp.*sin_x2;
            f4_6 = f4_6 .* 12 .* beta;
            
            f4 = (f4_1 + f4_2 + f4_3 + f4_4 + f4_5 + f4_6) ./ (4 * beta .* gamma);
            
            % State derivative
            dxdt(1, :) = x3;
            dxdt(2, :) = x4;
            dxdt(3, :) = f3;
            dxdt(4, :) = f4;
        end
        
        function [Ad, Bd, cd] = localLinearization(obj, xbar, ubar)
            % localLinearization: Calculates exact numerical Jacobians using Complex-Step Differentiation
            % Best practice alternative to Symbolic Math Toolbox for nonlinear models
            
            h = 1e-8;
            Ac = zeros(obj.n, obj.n);
            Bc = zeros(obj.n, obj.m);
            
            % State Jacobian A
            for i = 1:obj.n
                x_c = xbar;
                x_c(i) = x_c(i) + 1i * h;
                Ac(:, i) = imag(obj.continuousDynamics(x_c, ubar)) / h;
            end
            
            % Input Jacobian B
            u_c = ubar + 1i * h;
            Bc(:, 1) = imag(obj.continuousDynamics(xbar, u_c)) / h;
            
            % Constant term c
            c_val = obj.continuousDynamics(xbar, ubar);
            
            % Exact discretization via matrix exponential
            M = [Ac, Bc, c_val; 
                 zeros(2, obj.n + obj.m + 1)];
            ABc_d = expm(M * obj.deltaT);
            
            Ad = ABc_d(1:obj.n, 1:obj.n);
            Bd = ABc_d(1:obj.n, obj.n+1);
            cd = ABc_d(1:obj.n, obj.n+2);
        end
    end
end
