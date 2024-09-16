classdef TerminalKoopmanMPC < handle
    % TERMINALKOOPMANMPC Nonlinear MPC using fmincon for data-driven surrogates.
    % Incorporates state constraints, tightened bounds, and terminal conditions.
    
    properties
        predictor   % Data-driven surrogate (e.g., kEDMDPredictor)
        Np          % Prediction horizon
        Q, R        % Stage cost weights
        umin, umax  % Input constraints
        S_bounds    % State bounds [xmin, xmax] as [n x 2] matrix
        
        % Robustness parameters
        eta         % Uniform error bound
        Lbar        % Lipschitz constant of surrogate model
        
        % Terminal conditions
        P           % Terminal cost matrix (from LQR on linearized model)
        terminal_c  % Terminal constraint level (x_N' P x_N <= c)
    end
    
    methods
        function obj = TerminalKoopmanMPC(predictor, Np, Q, R, ubounds, S_bounds)
            obj.predictor = predictor;
            obj.Np = Np;
            obj.Q = Q;
            obj.R = R;
            obj.umin = ubounds(1);
            obj.umax = ubounds(2);
            obj.S_bounds = S_bounds;
            
            % Default robustness and terminal conditions
            obj.eta = 0;
            obj.Lbar = 1;
            obj.P = Q;
            obj.terminal_c = inf;
        end
        
        function setRobustnessParameters(obj, eta, Lbar)
            obj.eta = eta;
            obj.Lbar = Lbar;
        end
        
        function setTerminalConditions(obj, P, terminal_c)
            obj.P = P;
            obj.terminal_c = terminal_c;
        end
        
        function [u_opt, exitflag, info] = solve(obj, x0)
            % Use fmincon to solve the nonlinear program via single-shooting
            
            n = length(x0);
            m = 1; % Single input assumed
            
            % Optimization variables: U = [u_0, u_1, ..., u_{Np-1}]
            U0 = zeros(m * obj.Np, 1);
            
            % Input bounds
            lb = repmat(obj.umin, m * obj.Np, 1);
            ub = repmat(obj.umax, m * obj.Np, 1);
            
            % Objective function
            costFunc = @(U) obj.computeCost(x0, U);
            
            % Nonlinear constraints function
            nonlconFunc = @(U) obj.computeConstraints(x0, U);
            
            % Optimization options
            options = optimoptions('fmincon', 'Display', 'off', ...
                'Algorithm', 'sqp', ...
                'SpecifyObjectiveGradient', false, ...
                'SpecifyConstraintGradient', false);
            
            [U_opt, ~, exitflag, output] = fmincon(costFunc, U0, [], [], [], [], lb, ub, nonlconFunc, options);
            
            if exitflag >= 1 || exitflag == -3
                % If feasible or within tolerance, use the first control action
                u_opt = U_opt(1);
            else
                % Fallback: output 0 if completely failed
                u_opt = 0;
            end
            
            info = output;
        end
        
        function J = computeCost(obj, x0, U)
            % Compute the stage and terminal cost
            x = x0;
            J = 0;
            
            for k = 1:obj.Np
                u = U(k);
                % Stage cost
                J = J + x' * obj.Q * x + u' * obj.R * u;
                
                % Simulate forward
                x = obj.predictor.predict(x, u);
            end
            
            % Terminal cost
            J = J + x' * obj.P * x;
        end
        
        function [c, ceq] = computeConstraints(obj, x0, U)
            % Computes inequality constraints c(U) <= 0
            % Computes equality constraints ceq(U) == 0
            
            n = length(x0);
            c = [];
            ceq = [];
            
            x = x0;
            
            for k = 1:obj.Np
                u = U(k);
                x = obj.predictor.predict(x, u);
                
                if k < obj.Np
                    % Tightened state constraints: x_k in S (-) B_{c_bar(k) * eta}
                    % S_bounds is [n x 2] matrix: [xmin, xmax]
                    % c_bar(k) = sum_{i=0}^{k-1} Lbar^i
                    if obj.Lbar == 1
                        c_bar = k;
                    else
                        c_bar = (1 - obj.Lbar^k) / (1 - obj.Lbar);
                    end
                    
                    tightening = c_bar * obj.eta;
                    
                    for i = 1:n
                        % x_k >= xmin + tightening  -->  xmin + tightening - x_k <= 0
                        c(end+1, 1) = obj.S_bounds(i, 1) + tightening - x(i);
                        
                        % x_k <= xmax - tightening  -->  x_k - xmax + tightening <= 0
                        c(end+1, 1) = x(i) - obj.S_bounds(i, 2) + tightening;
                    end
                end
            end
            
            % Terminal constraint: x_N' P x_N <= terminal_c
            % x_N' P x_N - terminal_c <= 0
            c(end+1, 1) = x' * obj.P * x - obj.terminal_c;
        end
    end
end
