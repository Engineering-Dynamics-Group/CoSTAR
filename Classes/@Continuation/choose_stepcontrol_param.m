% This function is a method of the Continuation class.
% It defines default values for the step_control_param property, which are used by step control.
%
%@obj: Continuation class object
function obj = choose_stepcontrol_param(obj)

if isempty(obj.step_control_param)                                      % If step_control_param was not defined by user: Default value is set

    switch obj.step_control

        case 'corrector_iterations'
            it_nom = 3;                                                 % Nominal number of corrector iterations
            alpha_max = 3/180*pi;                                       % Maximum permissible angle (used by additional constraint)
            % y_rel_max = 0.002;                                        % Maximum value of additional constraint option 1
            % pred_step_rel_max = 0.1;                                  % Maximum value of additional constraint option 2
            % If additional constraint option 1 or 2 shall be used, uncomment the relevant lines in step_control.m as well as in this function ...
            % and replace alpha_max by the relevant parameter in the line below
            obj.step_control_param = [it_nom, alpha_max];               % Setting the parameter array


        case 'norm_corrector'
            it_max = 3;                                                 % Maximal number of corrector iterations
            norm_corr_nom = 0.0025;                                     % Nominal value of norm "corrector-vector"
            obj.step_control_param = [it_max, norm_corr_nom];           % Setting the parameter array


        case 'angle'
            it_max = 3;                                                 % Maximal number of corrector iterations
            alpha_nom = 3/180*pi;                                       % Nominal angle between the direction vectors at the last two points
            obj.step_control_param = [it_max, alpha_nom];               % Setting the parameter array


        case 'combination'
            it_nom = 3;                                                 % Nominal number of corrector iterations
            alpha_nom = 3/180*pi;                                       % Nominal angle between the direction vectors at the last two points
            obj.step_control_param = [it_nom, alpha_nom];               % Setting the parameter array


        case 'pid'
            it_max = 3;                                                 % Nominal number of corrector iterations
            tol = 0.1;                                                  % Reference value
            % Valli, Elias, Carey, Coutinho (2009) - "PID adaptive control of incremental and arclength continuation": kP = 0.075; kI = 0.175; kD = 0.01;
            kP = 0.075;
            kI = 0.175;
            kD = 0.01;
            obj.step_control_param = [it_max, tol, kP, kI, kD];         % Setting the parameter array

    end

end

end