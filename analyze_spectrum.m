% ANALYZE_SPECTRUM.m
% Validates the theoretical properties of the learned Koopman operator.

% 1. Extract Eigenvalues and Left/Right Eigenvectors
[W, D] = eig(koop.Alift);
eigenvalues = diag(D);

% 2. Plot Discrete-Time Spectrum
figure('Name', 'Koopman Spectrum');
hold on; grid on; axis equal;
th = linspace(0, 2*pi, 100);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1.5); % Unit circle

% Plot EDMD eigenvalues
scatter(real(eigenvalues), imag(eigenvalues), 50, 'b', 'filled');
title('Eigenvalues of Lifted Operator A_{lift}');
xlabel('Real'); ylabel('Imaginary');
legend('Unit Circle', '\lambda (Koopman)');
% Insight: Stable limits cycles (like Van der Pol) should have dominant
% eigenvalues exactly ON the unit circle, with the rest strictly inside.

% 3. Extract and Visualize the Principal Koopman Eigenfunction
% Left eigenvectors 'V' correspond to the Koopman eigenfunctions: 
% phi(x) = v_i^T * z(x). Let's plot the slowest decaying eigenfunction.
[V_left, ~] = eig(koop.Alift');
[~, sort_idx] = sort(abs(eigenvalues), 'descend');
dominant_eig_vec = V_left(:, sort_idx(2)); % index 1 is usually the trivial eigenvalue 1

% Evaluate on a grid
[X1, X2] = meshgrid(linspace(-2, 2, 50), linspace(-2, 2, 50));
X_grid = [X1(:)'; X2(:)'];
Z_grid = koop.lift(X_grid);

% Eigenfunction evaluation
Phi_grid = dominant_eig_vec' * Z_grid;
Phi_mag = reshape(abs(Phi_grid), 50, 50);

figure;
surf(X1, X2, Phi_mag, 'EdgeColor', 'none');
title('Magnitude of Dominant Koopman Eigenfunction');
xlabel('x_1'); ylabel('x_2'); zlabel('|\phi(x)|');
colormap parula; colorbar;
% Insight: The level sets of this eigenfunction correspond to the isochrons 
% (invariant sets) of the nonlinear oscillator!