classdef KoopmanMPC < handle
    % KOOPMANMPC Formulates a dense Quadratic Program (QP) for fast 
    % optimal control using the learned linear Koopman operator.
    
    properties
        koop        % Trained KoopmanPredictor
        Np = 10     % Prediction Horizon
        Q, R        % Cost weights (State and Input)
        umin, umax  % Input constraints
        
        % Condensed QP Matrices
        H, Gamma, Phi, Qbar, Rbar
    end
    
    methods
        function obj = KoopmanMPC(koop, Np, Q, R, ubounds)
            obj.koop = koop;
            obj.Np = Np;
            obj.Q = Q;
            obj.R = R;
            obj.umin = ubounds(1);
            obj.umax = ubounds(2);
            obj.buildMatrices();
        end
        
        function buildMatrices(obj)
            % Maps original state cost Q to lifted state cost: z' * (C'*Q*C) * z
            Q_lift = obj.koop.Clift' * obj.Q * obj.koop.Clift;
            A = obj.koop.Alift;
            B = obj.koop.Blift;
            Nlift = obj.koop.lifting.Nlift;
            
            % Build dense Prediction Matrices: Z = Phi*z0 + Gamma*U
            obj.Phi = zeros(Nlift * obj.Np, Nlift);
            obj.Gamma = zeros(Nlift * obj.Np, obj.Np); % scalar input m=1
            
            for k = 1:obj.Np
                obj.Phi((k-1)*Nlift+1 : k*Nlift, :) = A^k;
                for j = 1:k
                    obj.Gamma((k-1)*Nlift+1 : k*Nlift, j) = A^(k-j) * B;
                end
            end
            
            % Block diagonal cost matrices
            obj.Qbar = kron(eye(obj.Np), Q_lift);
            obj.Rbar = kron(eye(obj.Np), obj.R);
            
            % QP Hessian: H = Gamma' * Qbar * Gamma + Rbar
            obj.H = obj.Gamma' * obj.Qbar * obj.Gamma + obj.Rbar;
            obj.H = (obj.H + obj.H') / 2; % Ensure strict symmetry for quadprog
        end
        
        function [u_opt, exitflag] = solve(obj, x0, xref)
            % Solves the QP: min 0.5 * U'*H*U + f'*U 
            % subject to Umin <= U <= Umax
            
            z0 = obj.koop.lift(x0);
            
            % Reference in lifted space (assuming constant setpoint)
            Zref = repmat(obj.koop.lift(xref), obj.Np, 1);
            
            % Linear cost term: f = 2 * Gamma' * Qbar * (Phi*z0 - Zref)
            % (Dropping the factor of 2 because quadprog solves 0.5*x'Hx + f'x)
            f = obj.Gamma' * obj.Qbar * (obj.Phi * z0 - Zref);
            
            lb = repmat(obj.umin, obj.Np, 1);
            ub = repmat(obj.umax, obj.Np, 1);
            
            options = optimoptions('quadprog', 'Display', 'off');
            [U, ~, exitflag] = quadprog(obj.H, f, [], [], [], [], lb, ub, [], options);
            
            if exitflag >= 1
                u_opt = U(1); % Apply first control move (Receding Horizon)
            else
                u_opt = 0; % Fallback
            end
        end
    end
end