classdef KoopmanPredictor < handle
    % KOOPMANPREDICTOR Data-driven EDMD pipeline for extracting Koopman operator.
    % Uses Thin-Plate Spline RBFs alongside the original states.
    
    properties
        model       % Instance of VanDerPolModel to gather data from
        Nrbf        % Number of Radial Basis Functions
        cent        % RBF centers
        
        Alift       % Lifted discrete A matrix
        Blift       % Lifted discrete B matrix
        Clift       % Lifted discrete C matrix (projects back to state space)
    end
    
    methods
        function obj = KoopmanPredictor(model, num_rbf)
            obj.model = model;
            obj.Nrbf = num_rbf;
            
            % Centers are chosen randomly with uniform distribution on the unit box [-1, 1]^n
            obj.cent = rand(model.n, obj.Nrbf) * 2 - 1;
        end
        
        function train(obj, Nsim, Ntraj)
            % Set random seed for reproducibility in training data
            rng(42);
            
            disp('Starting data collection...');
            % Random forcing in [-1, 1]
            Ubig = 2 * rand(Nsim, Ntraj) - 1;
            
            % Random initial conditions in [-1, 1]
            Xcurrent = 2 * rand(obj.model.n, Ntraj) - 1;
            
            % Preallocate for speed
            total_samples = Nsim * Ntraj;
            X = zeros(obj.model.n, total_samples);
            Y = zeros(obj.model.n, total_samples);
            U = zeros(obj.model.m, total_samples);
            
            idx = 1;
            for i = 1:Nsim
                U_i = Ubig(i, :);
                Xnext = obj.model.discreteStep(Xcurrent, U_i);
                
                range = idx:(idx + Ntraj - 1);
                X(:, range) = Xcurrent;
                Y(:, range) = Xnext;
                U(:, range) = U_i;
                
                Xcurrent = Xnext;
                idx = idx + Ntraj;
            end
            disp('Data collection DONE.');
            
            disp('Starting LIFTING...');
            Xlift = obj.lift(X);
            Ylift = obj.lift(Y);
            disp('Lifting DONE.');
            
            disp('Starting REGRESSION (EDMD)...');
            W = [Ylift; X];
            V = [Xlift; U];
            
            % Regularization to prevent ill-conditioning
            VVt = V * V';
            WVt = W * V';
            % pseudo-inverse for robust least squares
            M = WVt * pinv(VVt);
            
            % Extract Koopman operator and projection matrices
            Nlift = obj.model.n + obj.Nrbf;
            
            obj.Alift = M(1:Nlift, 1:Nlift);
            obj.Blift = M(1:Nlift, Nlift+1:end);
            obj.Clift = M(Nlift+1:end, 1:Nlift);
            disp('Regression DONE. Koopman Predictor is ready.');
        end
        
        function x_lift = lift(obj, x)
            % lift: Transforms state x into lifted state space [x; phi(x)]
            % x: [n x N] matrix
            N = size(x, 2);
            phi = zeros(obj.Nrbf, N);
            
            % Compute Thin-Plate Spline RBF for each center
            for i = 1:obj.Nrbf
                c = obj.cent(:, i);
                % Euclidean distance to center
                % norm across columns using vectorized operations
                r = sqrt(sum((x - c).^2, 1));
                
                % r^2 * log(r)
                % To handle log(0) which returns -Inf, use r(idx) > 0 check
                idx = r > 0;
                phi(i, idx) = (r(idx).^2) .* log(r(idx));
            end
            
            x_lift = [x; phi];
        end
        
        function x_next = predict(obj, x_lift, u)
            % predict: advances lifted state by one step using linear
            % koopman matrices mapping.
            x_next = obj.Alift * x_lift + obj.Blift * u;
        end
        
        function x_original = project(obj, x_lift)
            % project: retrieves the original state variables from the
            % lifted space via the projection matrix C
            x_original = obj.Clift * x_lift;
        end
    end
end
