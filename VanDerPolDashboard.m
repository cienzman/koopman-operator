classdef VanDerPolDashboard < handle
    % VANDERPOLDASHBOARD Interactive dashboard for comparing trajectory predictors.
    %
    %   Three prediction strategies are compared against the ground-truth
    %   nonlinear simulation of the Van der Pol oscillator:
    %
    %     1. True (Nonlinear RK4)          — the reference simulation
    %     2. Koopman Predictor             — global linear model in lifted space
    %     3. Local Linearization at x(0)  — first-order Taylor expansion at the initial state
    %     4. Local Linearization at (0,0) — first-order Taylor expansion at the origin
    %
    %   Use the sliders to change the initial state and constant control input,
    %   then observe how well each approximation tracks the true trajectory.

    properties (Access = private)
        fig

        % Axes
        ax_phase        % Phase-plane portrait  (x2 vs x1)
        ax_x1           % Time series of x1
        ax_x2           % Time series of x2

        % Slider controls
        sl_x1           % Initial condition x1(0)
        sl_x2           % Initial condition x2(0)
        sl_u            % Constant control input u

        % Visibility toggles
        chk_true
        chk_koopman
        chk_loc_x0
        chk_loc_orig

        % Line handles (struct-of-structs, one per method)
        lph             % phase-plane lines
        lt1             % x1 time-series lines
        lt2             % x2 time-series lines

        % RMSE label handles
        lbl_rmse_koop
        lbl_rmse_locx0
        lbl_rmse_loc0

        % System objects
        vdp
        koop

        % Simulation length
        Nsim
        Tmax = 3;

        % Color palette (consistent across all axes)
        COLORS
    end

    methods (Access = public)
        function obj = VanDerPolDashboard(koopman_predictor)
            obj.koop  = koopman_predictor;
            obj.vdp   = koopman_predictor.model;
            obj.Nsim  = round(obj.Tmax / obj.vdp.deltaT);

            % Colorblind-friendly palette
            obj.COLORS = [0.00 0.45 0.70;   % blue   — True
                          0.90 0.60 0.00;   % orange — Koopman
                          0.00 0.62 0.45;   % green  — Loc x0
                          0.80 0.40 0.00];  % red    — Loc origin

            obj.buildUI();
            obj.runAndPlot();
        end
    end

    methods (Access = private)

        function buildUI(obj)
            obj.fig = uifigure( ...
                'Name',     'Koopman Predictor — Van der Pol Oscillator', ...
                'Position', [60 60 1340 700]);

            % ── Top-level layout: [plots | controls] ─────────────────────────
            root = uigridlayout(obj.fig, [1 2]);
            root.ColumnWidth = {'3x', '1x'};
            root.Padding     = [10 10 10 10];

            % ── Left: three stacked axes ──────────────────────────────────────
            left = uigridlayout(root, [3 1]);
            left.RowHeight  = {'1x', '1x', '1x'};
            left.RowSpacing = 8;

            obj.ax_phase = obj.makeAxes(left, ...
                'Phase Plane', 'x_1  (position)', 'x_2  (velocity)');
            obj.ax_x1 = obj.makeAxes(left, ...
                'State x_1 over time  (position)', 'Time  [s]', 'x_1');
            obj.ax_x2 = obj.makeAxes(left, ...
                'State x_2 over time  (velocity)', 'Time  [s]', 'x_2');

            % ── Right: control panel ──────────────────────────────────────────
            ctrl = uipanel(root, ...
                'Title',      'Controls', ...
                'FontSize',   12, ...
                'FontWeight', 'bold');

            cg = uigridlayout(ctrl, [14 1]);
            cg.RowHeight  = {22, 40, 22, 40, 22, 40, 14, ...
                             22, 22, 22, 22, 14, 80, 36};
            cg.Padding    = [12 10 12 10];
            cg.RowSpacing = 4;

            % Sliders
            obj.sl_x1 = obj.makeSlider(cg, 'Initial position  x1(0)', -2, 2,  0.5);
            obj.sl_x2 = obj.makeSlider(cg, 'Initial velocity  x2(0)', -2, 2,  0.5);
            obj.sl_u  = obj.makeSlider(cg, 'Control input  u  (constant)', -2, 2, 0.0);

            % Divider
            obj.makeDivider(cg);

            % Visibility checkboxes — colour-coded to match their lines
            obj.chk_true     = obj.makeCheckbox(cg, 'True nonlinear dynamics',  1, obj.COLORS(1,:));
            obj.chk_koopman  = obj.makeCheckbox(cg, 'Koopman predictor',        1, obj.COLORS(2,:));
            obj.chk_loc_x0   = obj.makeCheckbox(cg, 'Local linear  (at x0)',    1, obj.COLORS(3,:));
            obj.chk_loc_orig = obj.makeCheckbox(cg, 'Local linear  (at origin)',1, obj.COLORS(4,:));

            % Divider
            obj.makeDivider(cg);

            % RMSE readout panel
            rmse_pnl = uipanel(cg, 'Title', 'Prediction error (RMSE)', ...
                'FontSize', 10);
            rg = uigridlayout(rmse_pnl, [3 1]);
            rg.RowHeight  = {22, 22, 22};
            rg.Padding    = [6 4 6 4];
            rg.RowSpacing = 2;
            obj.lbl_rmse_koop  = uilabel(rg, 'Text', 'Koopman:      —', ...
                'FontColor', obj.COLORS(2,:), 'FontWeight', 'bold');
            obj.lbl_rmse_locx0 = uilabel(rg, 'Text', 'Local (x0):   —', ...
                'FontColor', obj.COLORS(3,:), 'FontWeight', 'bold');
            obj.lbl_rmse_loc0  = uilabel(rg, 'Text', 'Local (0,0):  —', ...
                'FontColor', obj.COLORS(4,:), 'FontWeight', 'bold');

            % Reset button
            uibutton(cg, 'Text', 'Reset to defaults', ...
                'ButtonPushedFcn', @(~,~) obj.reset(), ...
                'FontSize', 12);

            % ── Pre-create line objects on each axis ──────────────────────────
            keys   = {'true', 'koopman', 'loc_x0', 'loc_orig'};
            styles = {'-', '--', '-.', ':'};
            labels = {'True (nonlinear)', 'Koopman', 'Local lin. (x0)', 'Local lin. (0,0)'};

            for k = 1:4
                obj.lph.(keys{k}) = plot(obj.ax_phase, NaN, NaN, ...
                    styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 2, ...
                    'DisplayName', labels{k});
                obj.lt1.(keys{k}) = plot(obj.ax_x1, NaN, NaN, ...
                    styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 2, ...
                    'DisplayName', labels{k});
                obj.lt2.(keys{k}) = plot(obj.ax_x2, NaN, NaN, ...
                    styles{k}, 'Color', obj.COLORS(k,:), 'LineWidth', 2, ...
                    'DisplayName', labels{k});
            end

            legend(obj.ax_phase, 'Location', 'best', 'FontSize', 9);

            % ── Wire up callbacks after all controls exist ────────────────────
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
            x0    = [obj.sl_x1.Value; obj.sl_x2.Value];
            u_val = obj.sl_u.Value;

            N = obj.Nsim;
            X_true = zeros(2, N+1);  X_true(:,1) = x0;
            X_koop = zeros(2, N+1);  X_koop(:,1) = x0;
            X_lx0  = zeros(2, N+1);  X_lx0(:,1)  = x0;
            X_l0   = zeros(2, N+1);  X_l0(:,1)   = x0;

            % Lifted initial state for Koopman
            z = obj.koop.lift(x0);

            % Local linearizations computed once at t=0
            [Ad_x0, ~, cd_x0] = obj.vdp.localLinearization(x0,   u_val);
            [Ad_0,  ~, cd_0 ] = obj.vdp.localLinearization([0;0], u_val);

            for k = 1:N
                % True nonlinear RK4 step
                X_true(:, k+1) = obj.vdp.discreteStep(X_true(:, k), u_val);

                % Koopman: advance lifted state, then project back
                z              = obj.koop.predict(z, u_val);
                X_koop(:, k+1) = obj.koop.project(z);

                % Local linear at x0:  x+ = Ad*(x - x0) + cd + x0
                % Bd term is zero because u is held constant and equals ubar
                X_lx0(:, k+1) = Ad_x0 * (X_lx0(:, k) - x0) + cd_x0 + x0;

                % Local linear at origin:  x+ = Ad*x + cd
                X_l0(:, k+1)  = Ad_0  *  X_l0(:, k)          + cd_0;
            end

            t = (0:N) * obj.vdp.deltaT;

            % Phase plane
            obj.setXY(obj.lph.true,    X_true(1,:), X_true(2,:));
            obj.setXY(obj.lph.koopman, X_koop(1,:), X_koop(2,:));
            obj.setXY(obj.lph.loc_x0,  X_lx0(1,:),  X_lx0(2,:));
            obj.setXY(obj.lph.loc_orig,X_l0(1,:),   X_l0(2,:));

            % x1 time series
            obj.setXY(obj.lt1.true,    t, X_true(1,:));
            obj.setXY(obj.lt1.koopman, t, X_koop(1,:));
            obj.setXY(obj.lt1.loc_x0,  t, X_lx0(1,:));
            obj.setXY(obj.lt1.loc_orig,t, X_l0(1,:));

            % x2 time series
            obj.setXY(obj.lt2.true,    t, X_true(2,:));
            obj.setXY(obj.lt2.koopman, t, X_koop(2,:));
            obj.setXY(obj.lt2.loc_x0,  t, X_lx0(2,:));
            obj.setXY(obj.lt2.loc_orig,t, X_l0(2,:));

            % RMSE relative to true trajectory
            rmse = @(X) sqrt(mean(sum((X - X_true).^2, 1)));
            obj.lbl_rmse_koop.Text  = sprintf('Koopman:      %.4f', rmse(X_koop));
            obj.lbl_rmse_locx0.Text = sprintf('Local (x0):   %.4f', rmse(X_lx0));
            obj.lbl_rmse_loc0.Text  = sprintf('Local (0,0):  %.4f', rmse(X_l0));

            obj.applyVisibility();
        end

        % ── Visibility ────────────────────────────────────────────────────────

        function applyVisibility(obj)
            pairs = {obj.chk_true,     'true';
                     obj.chk_koopman,  'koopman';
                     obj.chk_loc_x0,   'loc_x0';
                     obj.chk_loc_orig, 'loc_orig'};
            for i = 1:size(pairs, 1)
                vis = matlab.lang.OnOffSwitchState(pairs{i,1}.Value);
                obj.lph.(pairs{i,2}).Visible = vis;
                obj.lt1.(pairs{i,2}).Visible = vis;
                obj.lt2.(pairs{i,2}).Visible = vis;
            end
        end

        % ── Reset ─────────────────────────────────────────────────────────────

        function reset(obj)
            obj.sl_x1.Value = 0.5;
            obj.sl_x2.Value = 0.5;
            obj.sl_u.Value  = 0.0;
            obj.runAndPlot();
        end

        % ── UI factory helpers ────────────────────────────────────────────────

        function ax = makeAxes(~, parent, ttl, xlbl, ylbl)
            ax = uiaxes(parent);
            title(ax,  ttl,  'FontSize', 11);
            xlabel(ax, xlbl, 'FontSize', 10);
            ylabel(ax, ylbl, 'FontSize', 10);
            grid(ax, 'on');
            hold(ax, 'on');
            ax.Box = 'on';
        end

        function sl = makeSlider(~, parent, label, lo, hi, val)
            uilabel(parent, 'Text', label, 'FontSize', 10);
            sl = uislider(parent, ...
                'Limits',      [lo hi], ...
                'Value',       val, ...
                'MajorTicks',  lo : (hi-lo)/4 : hi, ...
                'MinorTicks',  lo : (hi-lo)/20 : hi);
        end

        function chk = makeCheckbox(~, parent, label, val, color)
            chk = uicheckbox(parent, ...
                'Text',       label, ...
                'Value',      val, ...
                'FontColor',  color, ...
                'FontWeight', 'bold', ...
                'FontSize',   10);
        end

        function makeDivider(~, parent)
            uilabel(parent, 'Text', '', 'BackgroundColor', [0.8 0.8 0.8]);
        end

        function setXY(~, line, xd, yd)
            line.XData = xd;
            line.YData = yd;
        end

    end
end