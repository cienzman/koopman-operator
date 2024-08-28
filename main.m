% MAIN Script for Koopman Predictor Benchmarks
% Initializes the interactive multi-model dashboard

clear variables;
close all;
clc;

% Add all subfolders to path
addpath(genpath(pwd));

disp('========================================');
disp('      Koopman Operator Benchmarks       ');
disp('========================================');

disp('Initializing System and launching Dashboard...');
disp('This might take a few seconds to train the baseline Predictors...');

% The dashboard automatically trains Van der Pol, Duffing, and Inverted Pendulum models
app = KoopmanDashboard();

disp('Dashboard is ready! Select different systems from the dropdown menu.');