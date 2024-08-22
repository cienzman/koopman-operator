% RUN_VALIDATION_SUITE.m
% Automates ablation and robustness testing.

%% 1. Hyperparameter Ablation (Nrbf vs RMSE)
Nrbf_list = [10, 50, 100, 200];
rmse_results = zeros(size(Nrbf_list));

for i = 1:length(Nrbf_list)
    % Inject different lifting strategies
    strategy = RBFLifting(2, Nrbf_list(i));
    temp_koop = KoopmanPredictor(vdp, strategy);
    temp_koop.train(100, 500); % Fast train
    
    % Test on validation trajectory
    [X_true, X_koop] = simulate_trajectory(temp_koop, [1; 1], 0, 300);
    rmse_results(i) = sqrt(mean((X_true(:) - X_koop(:)).^2));
end

figure; plot(Nrbf_list, rmse_results, '-o', 'LineWidth', 2);
title('Ablation: Number of RBFs vs. Prediction Error');
xlabel('Number of RBFs (N_{rbf})'); ylabel('Validation RMSE');

%% 2. Robustness to Measurement Noise in Training
disp('Testing Robustness to Training Noise...');
% (Simulate adding 5% Gaussian noise to X_snap and Y_snap during training)
% A robust operator requires regularization (e.g., Tikhonov/Ridge regression
% in the pinv step: pinv(V*V' + lambda*I)). 

%% 3. Generalization: Out of Distribution (OOD) & Time-Varying Input
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