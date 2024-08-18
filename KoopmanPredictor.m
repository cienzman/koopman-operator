classdef KoopmanPredictor < handle
    % KOOPMANPREDICTOR Data-driven Koopman operator via EDMD.
    %
    %   The nonlinear state x is "lifted" into a higher-dimensional feature
    %   space using the original states plus Nrbf thin-plate-spline RBFs:
    %
    %       z = [x; phi(x)]     (lifted state, dimension n + Nrbf)
    %
    %   EDMD then fits linear matrices (Alift, Blift) such that
    %
    %       z_{k+1} ≈ Alift * z_k + Blift * u_k
    %
    %   and Clift recovers the original states:  x_k ≈ Clift * z_k.

    properties
        model       % VanDerPolModel instance
        Nrbf        % Number of RBF features
        cent        % RBF centers  [n x Nrbf]

        Alift       % Lifted state-transition matrix  [(n+Nrbf) x (n+Nrbf)]
        Blift       % Lifted input matrix             [(n+Nrbf) x m]
        Clift       % Projection back to state space  [n x (n+Nrbf)]
    end

    methods
        function obj = KoopmanPredictor(model, num_rbf)
            obj.model = model;
            obj.Nrbf  = num_rbf;
            % Distribute RBF centers uniformly in the state-space box [-1, 1]^n
            obj.cent  = rand(model.n, num_rbf) * 2 - 1;
        end

        function train(obj, Nsim, Ntraj)
            % Collect data, build lifted snapshot matrices, solve EDMD regression.
            rng(42);    % Fix seed so results are reproducible

            fprintf('Collecting %d training samples...\n', Nsim * Ntraj);
            U       = 2 * rand(Nsim, Ntraj) - 1;         % random inputs  in [-1, 1]
            Xcurrent = 2 * rand(obj.model.n, Ntraj) - 1; % random ICs     in [-1, 1]^n

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

            fprintf('Lifting snapshots...\n');
            Xlift = obj.lift(X_snap);
            Ylift = obj.lift(Y_snap);

            % EDMD least-squares:  solve  W ≈ M * V  for M
            %   rows of W: [lifted next-state; current state]
            %   rows of V: [lifted current-state; current input]
            fprintf('Solving EDMD regression...\n');
            W = [Ylift; X_snap];
            V = [Xlift; U_snap];
            M = W * V' * pinv(V * V');

            Nlift = obj.model.n + obj.Nrbf;
            obj.Alift = M(1:Nlift,        1:Nlift);
            obj.Blift = M(1:Nlift,        Nlift+1:end);
            obj.Clift = M(Nlift+1:end,    1:Nlift);
            fprintf('Training complete. Koopman predictor is ready.\n');
        end

        function z = lift(obj, x)
            % Maps state x -> lifted feature vector z = [x; phi(x)].
            % x : [n x N],   z : [(n + Nrbf) x N]
            N   = size(x, 2);
            phi = zeros(obj.Nrbf, N);

            for i = 1:obj.Nrbf
                % Euclidean distance from every sample to this RBF center
                r = sqrt(sum((x - obj.cent(:,i)).^2, 1));   % [1 x N]

                % Thin-plate-spline kernel:  r^2 * log(r)
                % Evaluated only where r > 0 to avoid log(0) = -Inf
                mask = r > 0;
                phi(i, mask) = r(mask).^2 .* log(r(mask));
            end

            z = [x; phi];
        end

        function z_next = predict(obj, z, u)
            % Advance the lifted state one step forward.
            z_next = obj.Alift * z + obj.Blift * u;
        end

        function x = project(obj, z)
            % Recover original state variables from the lifted state.
            x = obj.Clift * z;
        end
    end
end