% MAIN Script for Van der Pol Koopman Predictor
% This script initializes the system, trains the Koopman Model using EDMD,
% and opens the interactive dashboard.

clear variables;
close all;
clc;

disp('========================================');
disp('   Koopman Operator for Van der Pol     ');
disp('========================================');

% 1. Create Model
deltaT = 0.01;
vdp = VanDerPolModel(deltaT);

% 2. Initialize Koopman Predictor with 100 RBFs
num_rbf = 100;
koop = KoopmanPredictor(vdp, num_rbf);

% 3. Train Data
% Nsim = 200, Ntraj = 1000 => 200,000 samples
Nsim = 200;
Ntraj = 1000;
koop.train(Nsim, Ntraj);

% 4. Launch Interactive Dashboard
disp('Launching Dashboard...');
app = VanDerPolDashboard(koop);

disp('Dashboard is ready! Play with the sliders to change initial conditions.');
