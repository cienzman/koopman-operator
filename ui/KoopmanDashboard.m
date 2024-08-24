classdef KoopmanDashboard < handle
    % KOOPMANDASHBOARD Interactive visualization for Koopman Models.
    % Supports swapping between multiple dynamical models (e.g., Van der Pol, Duffing).
    
    properties
        fig
        ax_phase
        ax_time
        
        % Controls
        model_dropdown
        x1_slider
        x2_slider
        u_slider
        
        cb_true
        cb_koopman
        cb_loc_x0
        cb_loc_orig
        
        % Data
        active_model
        active_koop
        models_dict  % Structure containing initialized models and predictors
        
        % Simulation parameters
        Tmax = 3;
        Nsim = 300;
        
        % Plot lines
        lines_phase = struct();
        lines_time = struct();
    end
    
    methods
        function obj = KoopmanDashboard()
            obj.initModels();
            obj.buildUI();
            obj.updateSimulation();
        end
        
        function initModels(obj)
            % Initialize physical models
            deltaT = 0.01;
            vdp = VanDerPolModel(deltaT);
            duffing = DuffingModel(deltaT); % Default params
            
            % Initialize and train predictors
            disp('Training Van der Pol Predictor...');
            vdp_koop = KoopmanPredictor(vdp, RBFLifting(vdp.n, 100));
            vdp_koop.train(100, 500); % Fast train: 50,000 samples
            
            disp('Training Duffing Predictor...');
            duff_koop = KoopmanPredictor(duffing, RBFLifting(duffing.n, 100));
            duff_koop.train(100, 500);
            
            obj.models_dict.vdp.model = vdp;
            obj.models_dict.vdp.koop = vdp_koop;
            
            obj.models_dict.duffing.model = duffing;
            obj.models_dict.duffing.koop = duff_koop;
            
            % Set default
            obj.active_model = vdp;
            obj.active_koop = vdp_koop;
        end
        
        function switchModel(obj, model_name)
            if strcmp(model_name, 'Van der Pol Oscillator')
                obj.active_model = obj.models_dict.vdp.model;
                obj.active_koop = obj.models_dict.vdp.koop;
            elseif strcmp(model_name, 'Duffing Oscillator')
                obj.active_model = obj.models_dict.duffing.model;
                obj.active_koop = obj.models_dict.duffing.koop;
            end
            
            % Adjust simulation time if necessary or limits
            obj.updateSimulation();
        end
        
        function buildUI(obj)
            % Create main figure
            obj.fig = uifigure('Name', 'Koopman Operator Predictor Dashboard', ...
                               'Position', [100 100 1200 600]);
            
            gl = uigridlayout(obj.fig, [1, 2]);
            gl.ColumnWidth = {'1x', 300};
            
            % Left side: plots
            plot_gl = uigridlayout(gl, [2, 1]);
            obj.ax_phase = uiaxes(plot_gl);
            obj.ax_time = uiaxes(plot_gl);
            
            title(obj.ax_phase, 'Phase Plane: x_2 vs x_1');
            xlabel(obj.ax_phase, 'x_1');
            ylabel(obj.ax_phase, 'x_2');
            grid(obj.ax_phase, 'on'); hold(obj.ax_phase, 'on');
            
            title(obj.ax_time, 'Time Series: state over time');
            xlabel(obj.ax_time, 'Time [s]');
            ylabel(obj.ax_time, 'x_1, x_2');
            grid(obj.ax_time, 'on'); hold(obj.ax_time, 'on');
            
            % Right side: controls panel
            control_pnl = uipanel(gl, 'Title', 'Simulation Controls');
            c_gl = uigridlayout(control_pnl, [11, 1]);
            c_gl.RowHeight = {22, 22, 22, 22, 22, 22, 30, 22, 22, 22, 22};
            
            % Dropdown for model swapping
            uilabel(c_gl, 'Text', 'Dynamical System:');
            obj.model_dropdown = uidropdown(c_gl, ...
                'Items', {'Van der Pol Oscillator', 'Duffing Oscillator'}, ...
                'ValueChangedFcn', @(s,e) obj.switchModel(e.Value));
            
            % Sliders
            uilabel(c_gl, 'Text', 'Initial State: x_1(0)');
            obj.x1_slider = uislider(c_gl, 'Limits', [-2 2], 'Value', 0.5, ...
                'ValueChangedFcn', @(s,e) obj.updateSimulation());
                
            uilabel(c_gl, 'Text', 'Initial State: x_2(0)');
            obj.x2_slider = uislider(c_gl, 'Limits', [-2 2], 'Value', 0.5, ...
                'ValueChangedFcn', @(s,e) obj.updateSimulation());
                
            uilabel(c_gl, 'Text', 'Control Input: U');
            obj.u_slider = uislider(c_gl, 'Limits', [-2 2], 'Value', 0.0, ...
                'ValueChangedFcn', @(s,e) obj.updateSimulation());
                
            % Checkboxes
            obj.cb_true = uicheckbox(c_gl, 'Text', 'True Dynamics (Non-Linear)', 'Value', 1, ...
                'ValueChangedFcn', @(s,e) obj.refreshPlotsVisibility());
            obj.cb_koopman = uicheckbox(c_gl, 'Text', 'Koopman Predictor', 'Value', 1, ...
                'ValueChangedFcn', @(s,e) obj.refreshPlotsVisibility());
            obj.cb_loc_x0 = uicheckbox(c_gl, 'Text', 'Local Linearization (at x_0)', 'Value', 1, ...
                'ValueChangedFcn', @(s,e) obj.refreshPlotsVisibility());
            obj.cb_loc_orig = uicheckbox(c_gl, 'Text', 'Local Linearization (at Origin)', 'Value', 1, ...
                'ValueChangedFcn', @(s,e) obj.refreshPlotsVisibility());
                
            % Reset Button
            uibutton(c_gl, 'Text', 'Reset', 'ButtonPushedFcn', @(s,e) obj.resetControls());
            
            % Init plot lines
            colors = lines(4);
            obj.lines_phase.true    = plot(obj.ax_phase, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'True');
            obj.lines_phase.koopman = plot(obj.ax_phase, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Koopman');
            obj.lines_phase.loc_x0  = plot(obj.ax_phase, NaN, NaN, 'Color', colors(3,:), 'LineWidth', 2, 'LineStyle', '-.', 'DisplayName', 'Loc(x_0)');
            obj.lines_phase.loc_orig= plot(obj.ax_phase, NaN, NaN, 'Color', colors(4,:), 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', 'Loc(0)');
            
            obj.lines_time.true    = plot(obj.ax_time, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'True x_1');
            obj.lines_time.koopman = plot(obj.ax_time, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Koopman x_1');
            obj.lines_time.loc_x0  = plot(obj.ax_time, NaN, NaN, 'Color', colors(3,:), 'LineWidth', 2, 'LineStyle', '-.', 'DisplayName', 'Loc(x_0) x_1');
            obj.lines_time.loc_orig= plot(obj.ax_time, NaN, NaN, 'Color', colors(4,:), 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', 'Loc(0) x_1');
            
            legend(obj.ax_phase, 'Location', 'best');
            legend(obj.ax_time, 'Location', 'best');
        end
        
        function resetControls(obj)
            obj.x1_slider.Value = 0.5;
            obj.x2_slider.Value = 0.5;
            obj.u_slider.Value = 0.0;
            obj.updateSimulation();
        end
        
        function updateSimulation(obj)
            x0 = [obj.x1_slider.Value; obj.x2_slider.Value];
            u_val = obj.u_slider.Value;
            
            % Simulation arrays
            x_true = zeros(2, obj.Nsim + 1); x_true(:, 1) = x0;
            x_koop = zeros(2, obj.Nsim + 1); x_koop(:, 1) = x0;
            x_loc_x0 = zeros(2, obj.Nsim + 1); x_loc_x0(:, 1) = x0;
            x_loc_0 = zeros(2, obj.Nsim + 1); x_loc_0(:, 1) = x0;
            
            xlift = obj.active_koop.lift(x0);
            
            [Ad_x0, Bd_x0, cd_x0] = obj.active_model.localLinearization(x0, u_val);
            [Ad_0, Bd_0, cd_0] = obj.active_model.localLinearization([0;0], u_val);
            
            for k = 1:obj.Nsim
                x_true(:, k+1) = obj.active_model.discreteStep(x_true(:, k), u_val);
                
                xlift = obj.active_koop.predict(xlift, u_val);
                x_koop(:, k+1) = obj.active_koop.project(xlift);
                
                x_loc_x0(:, k+1) = Ad_x0 * (x_loc_x0(:, k) - x0) + Bd_x0 * (0) + cd_x0 + x0;
                x_loc_0(:, k+1) = Ad_0 * (x_loc_0(:, k) - 0) + Bd_0 * (0) + cd_0 + 0;
            end
            
            time_vec = (0:obj.Nsim) * obj.active_model.deltaT;
            
            obj.lines_phase.true.XData = x_true(1, :); obj.lines_phase.true.YData = x_true(2, :);
            obj.lines_phase.koopman.XData = x_koop(1, :); obj.lines_phase.koopman.YData = x_koop(2, :);
            obj.lines_phase.loc_x0.XData = x_loc_x0(1, :); obj.lines_phase.loc_x0.YData = x_loc_x0(2, :);
            obj.lines_phase.loc_orig.XData = x_loc_0(1, :); obj.lines_phase.loc_orig.YData = x_loc_0(2, :);
            
            obj.lines_time.true.XData = time_vec; obj.lines_time.true.YData = x_true(1, :);
            obj.lines_time.koopman.XData = time_vec; obj.lines_time.koopman.YData = x_koop(1, :);
            obj.lines_time.loc_x0.XData = time_vec; obj.lines_time.loc_x0.YData = x_loc_x0(1, :);
            obj.lines_time.loc_orig.XData = time_vec; obj.lines_time.loc_orig.YData = x_loc_0(1, :);
            
            axis(obj.ax_phase, 'auto');
            obj.refreshPlotsVisibility();
        end
        
        function refreshPlotsVisibility(obj)
            % Create a cell array of the two states, mapping 0->index 1, and 1->index 2
            state_opts = {'off', 'on'};
            on_off = @(val) state_opts{logical(val) + 1};
            
            obj.lines_phase.true.Visible = on_off(obj.cb_true.Value);
            obj.lines_time.true.Visible = on_off(obj.cb_true.Value);
            
            obj.lines_phase.koopman.Visible = on_off(obj.cb_koopman.Value);
            obj.lines_time.koopman.Visible = on_off(obj.cb_koopman.Value);
            
            obj.lines_phase.loc_x0.Visible = on_off(obj.cb_loc_x0.Value);
            obj.lines_time.loc_x0.Visible = on_off(obj.cb_loc_x0.Value);
            
            obj.lines_phase.loc_orig.Visible = on_off(obj.cb_loc_orig.Value);
            obj.lines_time.loc_orig.Visible = on_off(obj.cb_loc_orig.Value);
        end
    end
end