classdef Reaction_Kinetics_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        ComputationalMethodDropDown    matlab.ui.control.DropDown
        SelectmethodtosolvethereactionEquationLabel  matlab.ui.control.Label
        ReactionOrderSlider            matlab.ui.control.Slider
        ReactionOrderLabel             matlab.ui.control.Label
        PlotButton                     matlab.ui.control.Button
        InitialconcmolLEditField       matlab.ui.control.NumericEditField
        InitialconcmolLEditFieldLabel  matlab.ui.control.Label
        TotalTimeEditField             matlab.ui.control.NumericEditField
        TotalTimeEditFieldLabel        matlab.ui.control.Label
        TimeStepEditField              matlab.ui.control.NumericEditField
        TimeStepEditFieldLabel         matlab.ui.control.Label
        RateConstantKEditField         matlab.ui.control.NumericEditField
        RateConstantKEditFieldLabel    matlab.ui.control.Label
        UIAxes                         matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: PlotButton
        function PlotButtonPushed(app, event)
            k = app.RateConstantKEditField.Value;
    n = app.ReactionOrderSlider.Value;
    c = app.InitialconcmolLEditField.Value;
    t = app.TotalTimeEditField.Value;
    h = app.TimeStepEditField.Value;
    s = app.ComputationalMethodDropDown.Value;


    % To Check whether the value given by the user is evaluable or not
    if k <= 0
        uialert(app.UIFigure, 'Rate constant k must be > 0.', 'Invalid Input'); return;
    elseif n < 0
        uialert(app.UIFigure, 'Reaction order n must be non-negative.', 'Invalid Input'); return;
    elseif c <= 0
        uialert(app.UIFigure, 'Initial concentration must be > 0.', 'Invalid Input'); return;
    elseif t <= 0
        uialert(app.UIFigure, 'Final time must be > 0.', 'Invalid Input'); return;
    elseif h <= 0 || h >= t
        uialert(app.UIFigure, 'Time step h must be > 0 and < tf.', 'Invalid Input'); return;
    end

    % Function for dCa/dt
            function dCAdt = dCAdt(t,C,k,n)
                dCAdt = -k * C.^n;
            end
    t0 = 0;
    t = t0:h:t;
    N = length(t);
    Ca = zeros(1, N);
    Cb = zeros(1, N);
    Ca(1) = c;
    Cb(1) = 0;
                
if s=="Runge-Kutta 4"
    % RK4 Computation
    for i = 1:N-1
        k1 = h * dCAdt(t(i), Ca(i),k,n);
        k2 = h * dCAdt(t(i) + h/2, Ca(i) + h*(k1/2),k,n);
        k3 = h * dCAdt(t(i) + h/2, Ca(i) + h*(k2/2), k, n);
        k4 = h * dCAdt(t(i) + h, Ca(i) + k3, k , n);
        Ca(i+1) = Ca(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        if Ca(i+1) < 0
            Ca(i+1) = 0;
        end
        Cb(i+1) = c - Ca(i+1);
    end
end

if s == "Explicit euler"
    % Explicit Euler Computation
    for i = 1:N-1
        Ca(i+1) = Ca(i) + h * dCAdt(t(i), Ca(i), k, n);
        if Ca(i+1) < 0
            Ca(i+1) = 0;
        end
        Cb(i+1) = c - Ca(i+1);
    end
end

if s == "Implicit euler"
    % Implicit Euler Computation
    for i = 1:N-1
    % Initial guess
    Ca_guess = Ca(i);

    for iter = 1:20
        G = Ca_guess - Ca(i) + h * k * Ca_guess^n;
        dG = 1 + h * k * n * Ca_guess^(n - 1);

        % Newton-Raphson
        Ca_new = Ca_guess - G / dG;

        % Convergence check
        if abs(Ca_new - Ca_guess) < 1e-8
            break;
        end

        Ca_guess = Ca_new;
    end

    % Update the concentration
    Ca(i+1) = max(Ca_guess, 0);  % ensure non-negative
    Cb(i+1) = c - Ca(i+1);
    end
end


    cla(app.UIAxes);  % Clear axes
    plot(app.UIAxes, t, Ca, 'b-', 'LineWidth', 2); hold(app.UIAxes, 'on');
    plot(app.UIAxes, t, Cb, 'r--', 'LineWidth', 2);
    hold(app.UIAxes, 'off');

    % Labels and title
    legend(app.UIAxes, {'Reactant A', 'Product B'});
    grid(app.UIAxes, 'on');
    
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 334];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Reaction Kinetics')
            xlabel(app.UIAxes, 'Time')
            ylabel(app.UIAxes, 'Concentration')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [251 77 369 253];

            % Create RateConstantKEditFieldLabel
            app.RateConstantKEditFieldLabel = uilabel(app.UIFigure);
            app.RateConstantKEditFieldLabel.HorizontalAlignment = 'right';
            app.RateConstantKEditFieldLabel.Position = [18 294 93 22];
            app.RateConstantKEditFieldLabel.Text = 'Rate Constant K';

            % Create RateConstantKEditField
            app.RateConstantKEditField = uieditfield(app.UIFigure, 'numeric');
            app.RateConstantKEditField.Position = [126 294 100 22];

            % Create TimeStepEditFieldLabel
            app.TimeStepEditFieldLabel = uilabel(app.UIFigure);
            app.TimeStepEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeStepEditFieldLabel.Position = [24 213 60 22];
            app.TimeStepEditFieldLabel.Text = 'Time Step';

            % Create TimeStepEditField
            app.TimeStepEditField = uieditfield(app.UIFigure, 'numeric');
            app.TimeStepEditField.Position = [99 213 100 22];

            % Create TotalTimeEditFieldLabel
            app.TotalTimeEditFieldLabel = uilabel(app.UIFigure);
            app.TotalTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.TotalTimeEditFieldLabel.Position = [24 170 60 22];
            app.TotalTimeEditFieldLabel.Text = 'Total Time';

            % Create TotalTimeEditField
            app.TotalTimeEditField = uieditfield(app.UIFigure, 'numeric');
            app.TotalTimeEditField.Position = [99 170 100 22];

            % Create InitialconcmolLEditFieldLabel
            app.InitialconcmolLEditFieldLabel = uilabel(app.UIFigure);
            app.InitialconcmolLEditFieldLabel.HorizontalAlignment = 'right';
            app.InitialconcmolLEditFieldLabel.Position = [11 252 106 22];
            app.InitialconcmolLEditFieldLabel.Text = 'Initial conc. (mol/L)';

            % Create InitialconcmolLEditField
            app.InitialconcmolLEditField = uieditfield(app.UIFigure, 'numeric');
            app.InitialconcmolLEditField.Position = [132 252 100 22];

            % Create PlotButton
            app.PlotButton = uibutton(app.UIFigure, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed, true);
            app.PlotButton.Position = [62 77 100 23];
            app.PlotButton.Text = 'Plot';

            % Create ReactionOrderLabel
            app.ReactionOrderLabel = uilabel(app.UIFigure);
            app.ReactionOrderLabel.HorizontalAlignment = 'right';
            app.ReactionOrderLabel.Position = [90 29 86 22];
            app.ReactionOrderLabel.Text = 'Reaction Order';

            % Create ReactionOrderSlider
            app.ReactionOrderSlider = uislider(app.UIFigure);
            app.ReactionOrderSlider.Limits = [0 10];
            app.ReactionOrderSlider.Position = [239 45 325 3];

            % Create SelectmethodtosolvethereactionEquationLabel
            app.SelectmethodtosolvethereactionEquationLabel = uilabel(app.UIFigure);
            app.SelectmethodtosolvethereactionEquationLabel.HorizontalAlignment = 'center';
            app.SelectmethodtosolvethereactionEquationLabel.Position = [11 119 84 30];
            app.SelectmethodtosolvethereactionEquationLabel.Text = {'Computational'; 'Method'};

            % Create ComputationalMethodDropDown
            app.ComputationalMethodDropDown = uidropdown(app.UIFigure);
            app.ComputationalMethodDropDown.Items = {'Runge-Kutta 4', 'Explicit euler', 'Implicit euler'};
            app.ComputationalMethodDropDown.Position = [106 127 100 22];
            app.ComputationalMethodDropDown.Value = 'Runge-Kutta 4';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Reaction_Kinetics_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end