classdef VanDerPolDashboard < handle
    % VANDERPOLDASHBOARD Interactive dashboard for comparing trajectory predictors.
    % Upgraded with UX principles: Validity indicators, tooltips, error tracking,
    % and time-varying/MPC input profiles.

    properties (Access = private)
        fig
        
        % Axes
        ax_phase, ax_x1, ax_x2, ax_err
        
        % UI Controls
        sl_x1, sl_x2, sl_u
        dd_u_type       % Dropdown for input profile
        lamp_domain     % Lamp for Out-of-Distribution warning
        lbl_domain      % Label for the lamp
        
        % Visibility toggles
        chk_true, chk_koopman, chk_loc_x0, chk_loc_orig
        
        % Line handles (structs)
        lph, lt1, lt2, le
        
        % RMSE labels
        lbl_rmse_koop, lbl_rmse_locx0, lbl_rmse_loc0
        
        % System objects
        vdp, koop, mpc
        Nsim
        Tmax = 3;
        COLORS
    end

    methods (Access = public)
        function obj = VanDerPolDashboard(koopman_predictor)
            obj.koop  = koopman_predictor;
            obj.vdp   = koopman_predictor.model;
            obj.Nsim  = round(obj.Tmax / obj.vdp.deltaT);
            
            % Initialize MPC controller (Horizon=10, Q=diag([10, 10]), R=1)
            obj.mpc = KoopmanMPC(obj.koop, 10, eye(2)*10, 1, [-5, 5]);

            obj.COLORS = [0.00 0.45 0.70;   % blue   - True
                          0.90 0.60 0.00;   % orange - Koopman
                          0.00 0.62 0.45;   % green  - Loc x0
                          0.80 0.40 0.00];  % red    - Loc origin

            obj.buildUI();
            obj.runAndPlot();
        end
    end

    methods (Access = private)

        function buildUI(obj)
            % Increased height to 850 to comfortably fit 4 subplots
            obj.fig = uifigure('Name', 'Koopman Predictor Dashboard', 'Position', [60 60 1340 850]);

            root = uigridlayout(obj.fig, [1 2]);
            root.ColumnWidth = {'3x', '1x'};
            root.Padding     = [10 10 10 10];

            % ── Left: 4 stacked axes ──────────────────────────────────────────
            left = uigridlayout(root, [4 1]);
            left.RowHeight  = {'1x', '1x', '1x', '1x'};
            left.RowSpacing = 8;

            obj.ax_phase = obj.makeAxes(left, 'Phase Plane', 'x_1 (position)', 'x_2 (velocity)');
            obj.ax_x1    = obj.makeAxes(left, 'State x_1 (position)', 'Time [s]', 'x_1');
            obj.ax_x2    = obj.makeAxes(left, 'State x_2 (velocity)', 'Time [s]', 'x_2');
            obj.ax_err   = obj.makeAxes(left, 'Instantaneous Prediction Error ||x - x_{hat}||_2', 'Time [s]', 'Error Magnitude');

            % ── Right: control panel ──────────────────────────────────────────
            ctrl = uipanel(root, 'Title', 'Controls & Diagnostics', 'FontSize', 12, 'FontWeight', 'bold');

            cg = uigridlayout(ctrl, [17 1]);
            cg.RowHeight  = {22, 40, 22, 40, 22, 40, 22, 22, 14, 22, 22, 22, 22, 14, 80, 14, 36};
            cg.Padding    = [12 10 12 10];
            cg.RowSpacing = 4;

            % Initial Condition Sliders
            obj.sl_x1 = obj.makeSlider(cg, 'Initial position x1(0)', -2, 2, 0.5);
            obj.sl_x2 = obj.makeSlider(cg, 'Initial velocity x2(0)', -2, 2, 0.5);
            
            % Input Controls
            uilabel(cg, 'Text', 'Control Input Profile:', 'FontSize', 10, 'FontWeight', 'bold');
            obj.dd_u_type = uidropdown(cg, ...
                'Items', {'Constant', 'Step (t=1s)', 'Sine Wave (1 Hz)', 'Koopman MPC (Zero Ref)'}, ...
                'ValueChangedFcn', @(~,~) obj.runAndPlot());
            obj.sl_u = obj.makeSlider(cg, 'Input Amplitude / Constant u', -2, 2, 0.0);

            obj.makeDivider(cg);

            % Visibility Checkboxes + Educational Tooltips
            obj.chk_true     = obj.makeCheckbox(cg, 'True nonlinear dynamics',  1, obj.COLORS(1,:), 'Ground truth simulated via 4th-order Runge Kutta.');
            obj.chk_koopman  = obj.makeCheckbox(cg, 'Koopman predictor',        1, obj.COLORS(2,:), 'Global linear predictions in infinite-dimensional lifted feature space.');
            obj.chk_loc_x0   = obj.makeCheckbox(cg, 'Local linear (at x0)',     1, obj.COLORS(3,:), '1st-order Taylor Expansion around the starting point. Fails as state drifts.');
            obj.chk_loc_orig = obj.makeCheckbox(cg, 'Local linear (at origin)', 1, obj.COLORS(4,:), '1st-order Taylor Expansion around (0,0). Fails far from origin.');

            obj.makeDivider(cg);

            % RMSE readout panel
            rmse_pnl = uipanel(cg, 'Title', 'Prediction Error (Overall RMSE)', 'FontSize', 10);
            rg = uigridlayout(rmse_pnl, [3 1]);
            rg.RowHeight  = {22, 22, 22};
            rg.Padding    = [6 4 6 4];
            rg.RowSpacing = 2;
            obj.lbl_rmse_koop  = uilabel(rg, 'Text', 'Koopman:      —', 'FontColor', obj.COLORS(2,:), 'FontWeight', 'bold');
            obj.lbl_rmse_locx0 = uilabel(rg, 'Text', 'Local (x0):   —', 'FontColor', obj.COLORS(3,:), 'FontWeight', 'bold');
            obj.lbl_rmse_loc0  = uilabel(rg, 'Text', 'Local (0,0):  —', 'FontColor', obj.COLORS(4,:), 'FontWeight', 'bold');

            % Validity Lamp Layout
            lamp_grid = uigridlayout(cg, [1 2]);
            lamp_grid.ColumnWidth = {20, '1x'};
            lamp_grid.Padding = [0 0 0 0];
            obj.lamp_domain = uilamp(lamp_grid, 'Color', 'green');
            obj.lbl_domain  = uilabel(lamp_grid, 'Text', 'Valid: Inside Training Domain');

            % Reset button
            uibutton(cg, 'Text', 'Reset to defaults', 'ButtonPushedFcn', @(~,~) obj.reset(), 'FontSize', 12);

            % ── Pre-create line objects ───────────────────────────────────────
            keys   = {'true', 'koopman', 'loc_x0', 'loc_orig'};
            styles = {'-', '--', '-.', ':'};
            labels = {'True', 'Koopman', 'Local (x0)', 'Local (0,0)'};

            for k = 1:4
                obj.lph.(keys{k}) = plot(obj.ax_phase, NaN, NaN, styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 2, 'DisplayName', labels{k});
                obj.lt1.(keys{k}) = plot(obj.ax_x1, NaN, NaN, styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 2, 'DisplayName', labels{k});
                obj.lt2.(keys{k}) = plot(obj.ax_x2, NaN, NaN, styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 2, 'DisplayName', labels{k});
                
                % Error plots don't need 'True' line
                if k > 1
                    obj.le.(keys{k}) = plot(obj.ax_err, NaN, NaN, styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 1.5, 'DisplayName', labels{k});
                end
            end
            legend(obj.ax_phase, 'Location', 'best', 'FontSize', 9);
            legend(obj.ax_err, 'Location', 'northwest', 'FontSize', 9);

            % ── Wire Callbacks ────────────────────────────────────────────────
            obj.sl_x1.ValueChangedFcn        = @(~,~) obj.runAndPlot();
            obj.sl_x2.ValueChangedFcn        = @(~,~) obj.runAndPlot();
            obj.sl_u.ValueChangedFcn         = @(~,~) obj.runAndPlot();
            obj.chk_true.ValueChangedFcn     = @(~,~) obj.applyVisibility();
            obj.chk_koopman.ValueChangedFcn  = @(~,~) obj.applyVisibility();
            obj.chk_loc_x0.ValueChangedFcn   = @(~,~) obj.applyVisibility();
            obj.chk_loc_orig.ValueChangedFcn = @(~,~) obj.applyVisibility();
        end

        % ── Simulation ────────────────────────────────────────────────────────

        function runAndPlot(obj)
            x0 = [obj.sl_x1.Value; obj.sl_x2.Value];
            
            % Update Validity Lamp
            if abs(x0(1)) > 1.0 || abs(x0(2)) > 1.0
                obj.lamp_domain.Color = [0.8 0.1 0.1]; % Red
                obj.lbl_domain.Text = 'Warning: OOD! EDMD may degrade.';
            else
                obj.lamp_domain.Color = [0.2 0.8 0.2]; % Green
                obj.lbl_domain.Text = 'Valid: Inside Training Domain.';
            end

            N = obj.Nsim;
            X_true = zeros(2, N+1); X_true(:,1) = x0;
            X_koop = zeros(2, N+1); X_koop(:,1) = x0;
            X_lx0  = zeros(2, N+1); X_lx0(:,1)  = x0;
            X_l0   = zeros(2, N+1); X_l0(:,1)   = x0;
            
            U_seq = zeros(1, N); % Track input applied
            u_amp = obj.sl_u.Value;
            u_type = obj.dd_u_type.Value;

            z = obj.koop.lift(x0);

            % Linearize accurately including Bd term for time-varying inputs
            [Ad_x0, Bd_x0, cd_x0] = obj.vdp.localLinearization(x0, 0);
            [Ad_0,  Bd_0,  cd_0 ] = obj.vdp.localLinearization([0;0], 0);

            for k = 1:N
                t_k = (k-1) * obj.vdp.deltaT;
                
                % Determine input for this step based on UI selection
                if strcmp(u_type, 'Constant')
                    u_k = u_amp;
                elseif strcmp(u_type, 'Step (t=1s)')
                    u_k = u_amp * (t_k >= 1.0);
                elseif strcmp(u_type, 'Sine Wave (1 Hz)')
                    u_k = u_amp * sin(2 * pi * 1.0 * t_k);
                elseif strcmp(u_type, 'Koopman MPC (Zero Ref)')
                    u_k = obj.mpc.solve(X_true(:,k), [0;0]); % Regulate to origin
                end
                
                U_seq(k) = u_k;

                % True nonlinear RK4 step
                X_true(:, k+1) = obj.vdp.discreteStep(X_true(:, k), u_k);

                % Koopman EDMD Step
                z              = obj.koop.predict(z, u_k);
                X_koop(:, k+1) = obj.koop.project(z);

                % Local Linear Models
                X_lx0(:, k+1) = Ad_x0 * (X_lx0(:, k) - x0) + Bd_x0 * u_k + cd_x0 + x0;
                X_l0(:, k+1)  = Ad_0  *  X_l0(:, k)        + Bd_0  * u_k + cd_0;
            end

            t = (0:N) * obj.vdp.deltaT;

            % Update UI Plots
            obj.setXY(obj.lph.true,    X_true(1,:), X_true(2,:));
            obj.setXY(obj.lph.koopman, X_koop(1,:), X_koop(2,:));
            obj.setXY(obj.lph.loc_x0,  X_lx0(1,:),  X_lx0(2,:));
            obj.setXY(obj.lph.loc_orig,X_l0(1,:),   X_l0(2,:));

            obj.setXY(obj.lt1.true,    t, X_true(1,:));
            obj.setXY(obj.lt1.koopman, t, X_koop(1,:));
            obj.setXY(obj.lt1.loc_x0,  t, X_lx0(1,:));
            obj.setXY(obj.lt1.loc_orig,t, X_l0(1,:));

            obj.setXY(obj.lt2.true,    t, X_true(2,:));
            obj.setXY(obj.lt2.koopman, t, X_koop(2,:));
            obj.setXY(obj.lt2.loc_x0,  t, X_lx0(2,:));
            obj.setXY(obj.lt2.loc_orig,t, X_l0(2,:));
            
            % Compute and plot Instantaneous Errors ||x(t) - x_hat(t)||
            err_koop = sqrt(sum((X_koop - X_true).^2, 1));
            err_lx0  = sqrt(sum((X_lx0 - X_true).^2, 1));
            err_l0   = sqrt(sum((X_l0 - X_true).^2, 1));
            
            obj.setXY(obj.le.koopman,  t, err_koop);
            obj.setXY(obj.le.loc_x0,   t, err_lx0);
            obj.setXY(obj.le.loc_orig, t, err_l0);

            % Update Global RMSE Readouts
            obj.lbl_rmse_koop.Text  = sprintf('Koopman:      %.4f', mean(err_koop));
            obj.lbl_rmse_locx0.Text = sprintf('Local (x0):   %.4f', mean(err_lx0));
            obj.lbl_rmse_loc0.Text  = sprintf('Local (0,0):  %.4f', mean(err_l0));

            obj.applyVisibility();
        end

        function applyVisibility(obj)
            pairs = {obj.chk_true, 'true'; obj.chk_koopman, 'koopman';
                     obj.chk_loc_x0, 'loc_x0'; obj.chk_loc_orig, 'loc_orig'};
            for i = 1:size(pairs, 1)
                vis = matlab.lang.OnOffSwitchState(pairs{i,1}.Value);
                obj.lph.(pairs{i,2}).Visible = vis;
                obj.lt1.(pairs{i,2}).Visible = vis;
                obj.lt2.(pairs{i,2}).Visible = vis;
                
                if i > 1 % True model doesn't have an error plot line
                    obj.le.(pairs{i,2}).Visible = vis;
                end
            end
        end

        function reset(obj)
            obj.sl_x1.Value = 0.5;
            obj.sl_x2.Value = 0.5;
            obj.sl_u.Value  = 0.0;
            obj.dd_u_type.Value = 'Constant';
            obj.runAndPlot();
        end

        function ax = makeAxes(~, parent, ttl, xlbl, ylbl)
            ax = uiaxes(parent);
            title(ax, ttl, 'FontSize', 11);
            xlabel(ax, xlbl, 'FontSize', 10);
            ylabel(ax, ylbl, 'FontSize', 10);
            grid(ax, 'on'); hold(ax, 'on'); ax.Box = 'on';
        end

        function sl = makeSlider(~, parent, label, lo, hi, val)
            uilabel(parent, 'Text', label, 'FontSize', 10);
            sl = uislider(parent, 'Limits', [lo hi], 'Value', val, ...
                'MajorTicks', lo:(hi-lo)/4:hi, 'MinorTicks', lo:(hi-lo)/20:hi);
        end

        function chk = makeCheckbox(~, parent, label, val, color, tooltip)
            chk = uicheckbox(parent, 'Text', label, 'Value', val, ...
                'FontColor', color, 'FontWeight', 'bold', 'FontSize', 10, 'Tooltip', tooltip);
        end

        function makeDivider(~, parent)
            uilabel(parent, 'Text', '', 'BackgroundColor', [0.8 0.8 0.8]);
        end

        function setXY(~, line, xd, yd)
            line.XData = xd; line.YData = yd;
        end
    end
end