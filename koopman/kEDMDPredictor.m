classdef kEDMDPredictor < handle
    % KEDMDPREDICTOR Data-driven Kernel EDMD for Control Affine Systems
    % Implements the PI-kEDMD algorithm from "Data-driven MPC formulation"
    %
    % Assumes system structure: x+ = f(x, u) = g0(x) + G(x)u
    
    properties
        model       % Model instance
        d           % Number of virtual cluster points
        d_samples   % Samples per cluster
        
        X_centers   % Virtual observation points (d columns)
        g0_tilde    % Local drift approximations [n x d]
        G_tilde     % Local control approximations [n x m x d]
        
        K_inv       % Inverse kernel matrix [d x d]
        K0_hat      % Hat matrix for drift [d x d]
        K_hat       % Hat matrices for control [d x d x m]
        
        kernel_type % 'wendland' or 'tps'
        kernel_radius % Scaling radius for the Wendland kernel
        use_pi_constraint % boolean, defaults to true
    end
    
    methods
        function obj = kEDMDPredictor(model, d, d_samples, kernel_type)
            obj.model = model;
            obj.d = d;
            obj.d_samples = d_samples;
            if nargin > 3
                obj.kernel_type = kernel_type;
            else
                obj.kernel_type = 'tps'; % Defaulting to Thin-Plate Splines for better global support
            end
            obj.kernel_radius = 2.0; % Default Wendland support radius
            obj.use_pi_constraint = true;
        end
        
        function train(obj)
            rng(42); % Reproducibility
            
            n = obj.model.n;
            m = obj.model.m;
            
            fprintf('Training PI-kEDMD Predictor (%d clusters, %d samples each)...\n', obj.d, obj.d_samples);
            
            % 1. Create cluster centers. We must ensure x_1 = 0
            obj.X_centers = 2 * rand(n, obj.d) - 1;
            obj.X_centers(:, 1) = 0; % PI Constraint!
            
            obj.g0_tilde = zeros(n, obj.d);
            obj.G_tilde  = zeros(n, m, obj.d);
            
            % Heuristic local sampling radius
            r_X = 2 / obj.d^(1/n); 
            
            % Step 1: Data Preparation & Local Regression
            fprintf('Performing local regression at clusters...\n');
            for i = 1:obj.d
                x_i = obj.X_centers(:, i);
                
                % Sample local data triplets (x_ij, u_ij, x_next_ij)
                X_local = x_i + r_X * (2 * rand(n, obj.d_samples) - 1);
                U_local = 2 * rand(m, obj.d_samples) - 1;
                X_next_local = zeros(n, obj.d_samples);
                
                for j = 1:obj.d_samples
                    X_next_local(:, j) = obj.model.discreteStep(X_local(:, j), U_local(:, j));
                end
                
                if i == 1 && obj.use_pi_constraint
                    % Physics-Informed at origin: g0(0) = 0
                    obj.g0_tilde(:, i) = 0;
                    
                    % Reduced regression for G(0):
                    % min || [x_next] - G_1 * U_local ||_F
                    G_1 = X_next_local * pinv(U_local);
                    obj.G_tilde(:, :, i) = G_1;
                else
                    % Standard regression: min || [x_next] - [g0 G] * [1; U] ||_F
                    U_matrix = [ones(1, obj.d_samples); U_local];
                    H_i = X_next_local * pinv(U_matrix);
                    
                    obj.g0_tilde(:, i) = H_i(:, 1);
                    obj.G_tilde(:, :, i) = H_i(:, 2:end);
                end
            end
            
            % Step 2: Interpolation via Kernel matrices
            fprintf('Building global kernel matrices...\n');
            
            % Build K_X
            K_X = zeros(obj.d, obj.d);
            for i = 1:obj.d
                for j = 1:obj.d
                    K_X(i, j) = obj.kernel(obj.X_centers(:, i), obj.X_centers(:, j));
                end
            end
            
            % Add small regularization for numerical stability
            obj.K_inv = pinv(K_X + 1e-6 * eye(obj.d));
            
            % Build K0_hat
            K_g0 = zeros(obj.d, obj.d);
            for i = 1:obj.d
                for j = 1:obj.d
                    K_g0(i, j) = obj.kernel(obj.X_centers(:, i), obj.g0_tilde(:, j));
                end
            end
            obj.K0_hat = obj.K_inv * K_g0 * obj.K_inv;
            
            % Build K_hat for each control input
            obj.K_hat = zeros(obj.d, obj.d, m);
            for k = 1:m
                K_Gk = zeros(obj.d, obj.d);
                for i = 1:obj.d
                    for j = 1:obj.d
                        K_Gk(i, j) = obj.kernel(obj.X_centers(:, i), obj.G_tilde(:, k, j));
                    end
                end
                obj.K_hat(:, :, k) = obj.K_inv * K_Gk * obj.K_inv;
            end
            
            fprintf('PI-kEDMD Training complete.\n');
        end
        
        function val = kernel(obj, x, y)
            r = norm(x - y);
            if strcmp(obj.kernel_type, 'wendland')
                % Wendland C2 for 2D/3D (k=1)
                r_scaled = r / obj.kernel_radius;
                if r_scaled < 1
                    val = (1/20) * (1 - r_scaled)^4 * (4 * r_scaled + 1);
                else
                    val = 0;
                end
            else
                % Thin plate spline
                if r > 0
                    val = r^2 * log(r);
                else
                    val = 0;
                end
            end
        end
        
        function k_vec = get_k_vec(obj, x)
            k_vec = zeros(obj.d, 1);
            for i = 1:obj.d
                k_vec(i) = obj.kernel(obj.X_centers(:, i), x);
            end
        end
        
        function x_next = predict(obj, x, u)
            % Control-affine propagation: 
            % x+ = \sum_i [ (K0_hat x_X) + \sum_j u_j (K_j_hat x_X) ]_i Phi_{x_i}(x)
            
            n = obj.model.n;
            m = obj.model.m;
            x_next = zeros(n, 1);
            
            k_vec = obj.get_k_vec(x); % [d x 1] Evaluated kernels at x
            
            for dim = 1:n
                % Coordinate functions as observables: psi_l(x) = x_l
                psi_X = obj.X_centers(dim, :)'; % [d x 1] coordinate evaluations at cluster centers
                
                % Compute Koopman weights for this dimension
                weights = obj.K0_hat * psi_X;
                for j = 1:m
                    weights = weights + u(j) * (obj.K_hat(:,:,j) * psi_X);
                end
                
                % Inner product with kernel vector
                x_next(dim) = weights' * k_vec;
            end
        end
    end
end
