classdef KoopmanPredictor < handle
    % KOOPMANPREDICTOR Data-driven Koopman operator via EDMD.
    % Now utilizes the Strategy Pattern for flexible Lifting selection.
    
    properties
        model       % VanDerPolModel instance
        lifting     % Instance of LiftingStrategy (e.g., RBFLifting)
        
        Alift       % Lifted state-transition matrix  [(n+Nlift) x (n+Nlift)]
        Blift       % Lifted input matrix             [(n+Nlift) x m]
        Clift       % Projection back to state space  [n x (n+Nlift)]
    end
    
    methods
        function obj = KoopmanPredictor(model, lifting_strategy)
            obj.model = model;
            obj.lifting = lifting_strategy;
        end
        
        function train(obj, Nsim, Ntraj)
            rng(42); % Fix seed for reproducibility
            
            fprintf('Collecting %d training samples...\n', Nsim * Ntraj);
            U        = 2 * rand(Nsim, Ntraj) - 1;         % Random inputs in [-1, 1]
            Xcurrent = 2 * rand(obj.model.n, Ntraj) - 1;  % Random ICs in [-1, 1]^n
            
            % Pre-allocate snapshot matrices
            X_snap = zeros(obj.model.n, Nsim * Ntraj);
            Y_snap = zeros(obj.model.n, Nsim * Ntraj);
            U_snap = zeros(obj.model.m, Nsim * Ntraj);
            
            for i = 1:Nsim
                u_i   = U(i, :);
                Xnext = obj.model.discreteStep(Xcurrent, u_i);
                
                cols = (i-1)*Ntraj + (1:Ntraj);
                X_snap(:, cols) = Xcurrent;
                Y_snap(:, cols) = Xnext;
                U_snap(:, cols) = u_i;
                
                Xcurrent = Xnext;
            end
            fprintf('Data collection done.\n');
            
            % 1. Fit the lifting strategy (e.g., K-Means centers for RBFs)
            fprintf('Fitting lifting strategy...\n');
            obj.lifting.fit(X_snap);
            
            % 2. Lift Snapshots
            fprintf('Lifting snapshots...\n');
            Xlift = obj.lifting.lift(X_snap);
            Ylift = obj.lifting.lift(Y_snap);
            
            % 3. EDMD least-squares regression
            fprintf('Solving EDMD regression...\n');
            W = [Ylift; X_snap];
            V = [Xlift; U_snap];
            M = W * V' * pinv(V * V');
            
            % Extract A, B, C matrices
            Nlift_total = obj.lifting.Nlift;
            obj.Alift = M(1:Nlift_total,        1:Nlift_total);
            obj.Blift = M(1:Nlift_total,        Nlift_total+1:end);
            obj.Clift = M(Nlift_total+1:end,    1:Nlift_total);
            
            fprintf('Training complete. Koopman predictor is ready.\n');
        end
        
        % --- Forward to Lifting Strategy ---
        function z = lift(obj, x)
            z = obj.lifting.lift(x);
        end
        
        function z_next = predict(obj, z, u)
            z_next = obj.Alift * z + obj.Blift * u;
        end
        
        function x = project(obj, z)
            x = obj.Clift * z;
        end
    end
end