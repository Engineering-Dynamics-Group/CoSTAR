% This function is a method of the Continuation class.
% It defines default values for the step_control_param property, which are used by step control.
%
% @obj: Continuation class object

function obj = choose_stepcontrol_param(obj)

if isempty(obj.step_control_param)                                      % If step_control_param was not defined by user: Default value is set

    switch obj.step_control

        case 'corrector_iterations'
            it_nom = 3;                                                 % Nominal number of corrector iterations
            add_cstr_nom = 3;                                           % Nominal angle in degrees (used by additional constraint)
            % add_cstr_nom = 0.003;                                     % Nominal value of additional constraint option 1
            % If additional constraint option 1 shall be used, uncomment the relevant lines in step_control.m as well as in this function
            obj.step_control_param = [it_nom, add_cstr_nom];            % Setting the parameter array


        case 'norm_corrector'
            it_nom = 2;                                                 % Nominal number of corrector iterations
            norm_corr_nom = 0.0025;                                     % Nominal value of norm "corrector-vector"
            obj.step_control_param = [it_nom, norm_corr_nom];           % Setting the parameter array


        case 'angle'
            it_nom = 2;                                                 % Nominal number of corrector iterations
            alpha_nom = 3;                                              % Nominal angle in degrees between the direction vectors at the last two points
            obj.step_control_param = [it_nom, alpha_nom];               % Setting the parameter array


        case 'combination'
            it_nom = 2;                                                 % Nominal number of corrector iterations
            alpha_nom = 3;                                              % Nominal angle in degrees between the direction vectors at the last two points
            obj.step_control_param = [it_nom, alpha_nom];               % Setting the parameter array


        case 'pid'
            it_nom = 2;                                                 % Nominal number of corrector iterations
            tol = 0.1;                                                  % Reference value
            % Valli, Elias, Carey, Coutinho (2009) - "PID adaptive control of incremental and arc-length continuation": kP = 0.075; kI = 0.175; kD = 0.01;
            kP = 0.075;
            kI = 0.175;
            kD = 0.01;
            obj.step_control_param = [it_nom, tol, kP, kI, kD];         % Setting the parameter array

    end

end

end