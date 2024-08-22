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

% 2. Setup Lifting Strategy & Initialize Predictor
num_rbf = 100;
lifting_strategy = RBFLifting(vdp.n, num_rbf);
koop = KoopmanPredictor(vdp, lifting_strategy);

% 3. Train Data
% Nsim = 200, Ntraj = 1000 => 200,000 samples
Nsim = 200;
Ntraj = 1000;
koop.train(Nsim, Ntraj);

% 4. Launch Interactive Dashboard
disp('Launching Dashboard...');
app = VanDerPolDashboard(koop);

disp('Dashboard is ready!');
disp('- Try setting the Input Profile to "Koopman MPC"');
disp('- Try moving the sliders to see the Valid/OOD Lamp change');