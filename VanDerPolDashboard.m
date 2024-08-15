classdef VanDerPolDashboard < handle
    % VANDERPOLDASHBOARD Interactive visualization for the Koopman Model.
    
    properties
        fig
        ax_phase
        ax_time
        
        % Controls
        x1_slider
        x2_slider
        u_slider
        
        cb_true
        cb_koopman
        cb_loc_x0
        cb_loc_orig
        
        % Models
        vdp
        koop
        
        % Simulation parameters
        Tmax = 3;
        Nsim = 300;
        
        % Plot lines
        lines_phase = struct();
        lines_time = struct();
    end
    
    methods
        function obj = VanDerPolDashboard(koopman_predictor)
            obj.koop = koopman_predictor;
            obj.vdp = koopman_predictor.model;
            obj.Nsim = round(obj.Tmax / obj.vdp.deltaT);
            
            obj.buildUI();
            obj.updateSimulation();
        end
        
        function buildUI(obj)
            % Create main figure
            obj.fig = uifigure('Name', 'Koopman Operator Predictor Dashboard', ...
                               'Position', [100 100 1200 600]);
            
            % Use a grid layout
            gl = uigridlayout(obj.fig, [1, 2]);
            gl.ColumnWidth = {'1x', 300};
            
            % Left side: plots
            plot_gl = uigridlayout(gl, [2, 1]);
            obj.ax_phase = uiaxes(plot_gl);
            obj.ax_time = uiaxes(plot_gl);
            
            title(obj.ax_phase, 'Phase Plane: x_2 vs x_1');
            xlabel(obj.ax_phase, 'x_1');
            ylabel(obj.ax_phase, 'x_2');
            grid(obj.ax_phase, 'on');
            hold(obj.ax_phase, 'on');
            
            title(obj.ax_time, 'Time Series: state over time');
            xlabel(obj.ax_time, 'Time [s]');
            ylabel(obj.ax_time, 'x_1, x_2');
            grid(obj.ax_time, 'on');
            hold(obj.ax_time, 'on');
            
            % Right side: controls panel
            control_pnl = uipanel(gl, 'Title', 'Simulation Controls');
            c_gl = uigridlayout(control_pnl, [9, 1]);
            c_gl.RowHeight = {22, 22, 22, 22, 30, 22, 22, 22, 22};
            
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
            
            % Phase plane lines
            obj.lines_phase.true    = plot(obj.ax_phase, NaN, NaN, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'True');
            obj.lines_phase.koopman = plot(obj.ax_phase, NaN, NaN, 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Koopman');
            obj.lines_phase.loc_x0  = plot(obj.ax_phase, NaN, NaN, 'Color', colors(3,:), 'LineWidth', 2, 'LineStyle', '-.', 'DisplayName', 'Loc(x_0)');
            obj.lines_phase.loc_orig= plot(obj.ax_phase, NaN, NaN, 'Color', colors(4,:), 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', 'Loc(0)');
            
            % Time series lines (just plotting x1 for simplicity, solid and dashed)
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
            
            % 1. True Dynamics
            x_true = zeros(2, obj.Nsim + 1);
            x_true(:, 1) = x0;
            
            % 2. Koopman Predictor
            xlift = obj.koop.lift(x0);
            x_koop = zeros(2, obj.Nsim + 1);
            x_koop(:, 1) = x0;
            
            % 3. Local Linearization at x0
            [Ad_x0, Bd_x0, cd_x0] = obj.vdp.localLinearization(x0, u_val);
            x_loc_x0 = zeros(2, obj.Nsim + 1);
            x_loc_x0(:, 1) = x0;
            
            % 4. Local Linearization at origin (0,0)
            [Ad_0, Bd_0, cd_0] = obj.vdp.localLinearization([0;0], u_val);
            x_loc_0 = zeros(2, obj.Nsim + 1);
            x_loc_0(:, 1) = x0;
            
            % Run Simulation Loop
            for k = 1:obj.Nsim
                % True
                x_true(:, k+1) = obj.vdp.discreteStep(x_true(:, k), u_val);
                
                % Koopman
                xlift = obj.koop.predict(xlift, u_val);
                x_koop(:, k+1) = obj.koop.project(xlift);
                
                % Loc x0  (x_{k+1} = Ad*(x_k - xbar) + Bd*(u_k - ubar) + cd + xbar)
                x_loc_x0(:, k+1) = Ad_x0 * (x_loc_x0(:, k) - x0) + Bd_x0 * (u_val - u_val) + cd_x0 + x0;
                
                % Loc 0
                x_loc_0(:, k+1) = Ad_0 * (x_loc_0(:, k) - 0) + Bd_0 * (u_val - u_val) + cd_0 + 0;
            end
            
            % Update plot data
            time_vec = (0:obj.Nsim) * obj.vdp.deltaT;
            
            % Phase plane
            obj.lines_phase.true.XData = x_true(1, :);
            obj.lines_phase.true.YData = x_true(2, :);
            
            obj.lines_phase.koopman.XData = x_koop(1, :);
            obj.lines_phase.koopman.YData = x_koop(2, :);
            
            obj.lines_phase.loc_x0.XData = x_loc_x0(1, :);
            obj.lines_phase.loc_x0.YData = x_loc_x0(2, :);
            
            obj.lines_phase.loc_orig.XData = x_loc_0(1, :);
            obj.lines_phase.loc_orig.YData = x_loc_0(2, :);
            
            % Time series
            obj.lines_time.true.XData = time_vec;
            obj.lines_time.true.YData = x_true(1, :);
            
            obj.lines_time.koopman.XData = time_vec;
            obj.lines_time.koopman.YData = x_koop(1, :);
            
            obj.lines_time.loc_x0.XData = time_vec;
            obj.lines_time.loc_x0.YData = x_loc_x0(1, :);
            
            obj.lines_time.loc_orig.XData = time_vec;
            obj.lines_time.loc_orig.YData = x_loc_0(1, :);
            
            % Auto-adjust limits slightly
            axis(obj.ax_phase, 'auto');
            
            obj.refreshPlotsVisibility();
        end
        
        function refreshPlotsVisibility(obj)
            % True
            if obj.cb_true.Value
                obj.lines_phase.true.Visible = 'on';
                obj.lines_time.true.Visible = 'on';
            else
                obj.lines_phase.true.Visible = 'off';
                obj.lines_time.true.Visible = 'off';
            end
            
            % Koopman
            if obj.cb_koopman.Value
                obj.lines_phase.koopman.Visible = 'on';
                obj.lines_time.koopman.Visible = 'on';
            else
                obj.lines_phase.koopman.Visible = 'off';
                obj.lines_time.koopman.Visible = 'off';
            end
            
            % Local x0
            if obj.cb_loc_x0.Value
                obj.lines_phase.loc_x0.Visible = 'on';
                obj.lines_time.loc_x0.Visible = 'on';
            else
                obj.lines_phase.loc_x0.Visible = 'off';
                obj.lines_time.loc_x0.Visible = 'off';
            end
            
            % Local 0
            if obj.cb_loc_orig.Value
                obj.lines_phase.loc_orig.Visible = 'on';
                obj.lines_time.loc_orig.Visible = 'on';
            else
                obj.lines_phase.loc_orig.Visible = 'off';
                obj.lines_time.loc_orig.Visible = 'off';
            end
        end
    end
end
