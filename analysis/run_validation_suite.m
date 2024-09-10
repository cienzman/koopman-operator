% RUN_VALIDATION_SUITE.m
% Automates ablation and robustness testing.

%% 1. Initialization
deltaT = 0.01;
vdp = VanDerPolModel(deltaT);

% Train a baseline Koopman Operator (100 RBFs)
disp('Training Baseline Koopman (100 RBFs)...');
baseline_strategy = RBFLifting(vdp.n, 100);
koop = KoopmanPredictor(vdp, baseline_strategy);
koop.train(200, 1000); % Baseline 200,000 samples

%% 2. Physics-Informed Proportional Error Bounds (Origin Equilibrium Test)
disp('Testing Physics-Informed Koopman at the Origin...');
pi_strategy = PIRBFLifting(vdp.n, 100);
koop_pi = KoopmanPredictor(vdp, pi_strategy);
koop_pi.train(200, 1000);

% Test prediction at exactly the origin
x0_origin = [0; 0];
u0 = 0;

z_std = koop.lift(x0_origin);
x_next_std = koop.project(koop.predict(z_std, u0));

z_pi = koop_pi.lift(x0_origin);
x_next_pi = koop_pi.project(koop_pi.predict(z_pi, u0));

fprintf('--- Equilibrium Test (x=0, u=0) ---\n');
fprintf('True Nonlinear model:      [%f; %f]\n', 0, 0);
fprintf('Standard Koopman Predictor:[%f; %f] (Error: %e)\n', x_next_std(1), x_next_std(2), norm(x_next_std));
fprintf('Physics-Informed Koopman:  [%f; %f] (Error: %e)\n', x_next_pi(1), x_next_pi(2), norm(x_next_pi));
fprintf('-----------------------------------\n\n');

%% 3. Hyperparameter Ablation (Nrbf vs RMSE)
disp('Running RBF Ablation Study...');
Nrbf_list = [10, 50, 100, 200];
rmse_results = zeros(size(Nrbf_list));

for i = 1:length(Nrbf_list)
    % Inject different lifting strategies
    strategy = RBFLifting(vdp.n, Nrbf_list(i));
    temp_koop = KoopmanPredictor(vdp, strategy);
    temp_koop.train(100, 500); % Fast train for ablation
    
    % Test on validation trajectory
    [X_true, X_koop] = simulate_trajectory(temp_koop, vdp, [1; 1], 0, 300);
    rmse_results(i) = sqrt(mean((X_true(:) - X_koop(:)).^2));
end

figure; plot(Nrbf_list, rmse_results, '-o', 'LineWidth', 2);
title('Ablation: Number of RBFs vs. Prediction Error');
xlabel('Number of RBFs (N_{rbf})'); ylabel('Validation RMSE');

%% 4. Robustness to Measurement Noise in Training
disp('Testing Robustness to Training Noise...');
% (Simulate adding 5% Gaussian noise to X_snap and Y_snap during training)
% A robust operator requires regularization (e.g., Tikhonov/Ridge regression
% in the pinv step: pinv(V*V' + lambda*I)). 

%% 5. Generalization: Out of Distribution (OOD) & Time-Varying Input
disp('Testing Time-Varying Input Tracking...');
Nsim = 500;
U_sin = sin(3 * (0:Nsim-1) * vdp.deltaT); % Sine wave input

x_curr_true = [1.5; -1.5]; % Starting outside [-1,1]
z_curr_koop = koop.lift(x_curr_true);

X_true = zeros(2, Nsim); X_koop = zeros(2, Nsim);

for k = 1:Nsim
    % True
    X_true(:,k) = x_curr_true;
    x_curr_true = vdp.discreteStep(x_curr_true, U_sin(k));
    
    % Koopman
    X_koop(:,k) = koop.project(z_curr_koop);
    z_curr_koop = koop.predict(z_curr_koop, U_sin(k));
end

figure;
plot(1:Nsim, X_true(1,:), 'b-', 'LineWidth', 2); hold on;
plot(1:Nsim, X_koop(1,:), 'r--', 'LineWidth', 2);
title('Generalization: Tracking Sine-Wave Input (Out of Distribution Initial State)');
legend('True Nonlinear', 'Koopman Predictor');

%% Helper Functions
function [X_true, X_koop] = simulate_trajectory(koop, model, x0, u0, Nsim)
    % Simulates a trajectory for Nsim steps under constant input u0
    X_true = zeros(model.n, Nsim);
    X_koop = zeros(model.n, Nsim);
    
    x_curr = x0;
    z_curr = koop.lift(x0);
    
    for k = 1:Nsim
        X_true(:, k) = x_curr;
        X_koop(:, k) = koop.project(z_curr);
        
        x_curr = model.discreteStep(x_curr, u0);
        z_curr = koop.predict(z_curr, u0);
    end
end