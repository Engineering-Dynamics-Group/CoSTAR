% Function step_control adapts step_width
%
% @obj: Continuation class object

function obj = stepcontrol(obj)

switch obj.step_control
 
    case 'off'                                  
        % No step control. Step width is a constant value (exception: corrector did not converge)
        % If step width was reduced because corrector did not converge: step width must be reset to initial step width
        % If step control is on: Step width does not have to be reset to initial step width, because step width is controlled via step control

        if obj.p_convergence == 0                                                                   % If corrector did not converge previously
            obj.step_width = obj.p_step_width_init;                                                 % Reset step-width to initial step width
            disp(append('Step width reset to step_width = ',num2str(obj.step_width)));              % Display information
        end


    otherwise
        
        r_limit = [0.5, 1.5];                                           % Limits of r. r is factor to adapt step width
        step_width_old = obj.step_width;                                % Store preceding step width (can be used to control display output)

        switch  obj.step_control

            case 'corrector_iterations'         
                % Step control based on number of corrector iterations

                % Step control parameters
                it_nom = obj.step_control_param(1);                     % Nominal number of corrector iterations
                add_cstr_value = obj.step_control_param(2);             % Parameter used by the additional constraint
        
                % Additional constraint: predictor point yp has to be near new point on curve p_y1
                % Option 1: ratio of norm of predictor point to norm of latest point on curve (activate the next two lines and deactivate option 3 to use this option)
                % y_rel = abs( norm(obj.yp) / norm(obj.y0) - 1);        % -1 in order to use a simple if-condition
                % if y_rel > add_cstr_value                             % If y_rel is too large: predictor point yp is too far away from y0
                % Option 2: ratio of norm "predictor vector" to norm of latest point on curve (activate the next two lines and deactivate option 3 to use this option). This can lead to step widths being more regular
                % pred_step_rel = obj.step_width / norm(obj.y0);        % norm(obj.step_width*obj.p_dy_old) = obj.step_width, as norm(obj.p_dy_old) = 1
                % if pred_step_rel > add_cstr_value                     % If pred_step_rel is too large: ???
                % Option 3: angle between direction vectors at the last two points on curve (this option is used by default)
                alpha = acos( obj.dy0' * obj.p_dy_old );                % Calculate angle between the last two direction vectors
                if alpha > add_cstr_value                               % If the angle is too large: discretization is too coarse with respect to local curvature
                    r_pre = r_limit(1);                                 % Minimal factor of r is taken to reduce step width
                    display_message_info = append('Angle between the last two direction vectors is too large: alpha = ', num2str(alpha), ' > ', num2str(add_cstr_value), ' = alpha_max. corr_it = ', num2str(obj.p_it));
        
                % Main criterion: adapt step width based on ratio of nominal number of corrector iterations to actual number of corrector iterations
                else
                    r_pre = it_nom / obj.p_it;                          % Calculate preliminary value of r using the ratio of nominal corrector iterations to actual corrector iterations
                    display_message_info = append('corr_it = ', num2str(obj.p_it));     % Set display message
                end
                d = 0.5;                                                % Set damping factor (default: 0.5)

            
            case 'norm_corrector'
                % Step control originally based on norm  of "corrector-vector" (= obj.p_y1 - obj.yp) to norm "predictor-vector" (= obj.dy0 * obj.step_width)
                % If the ratio is large, the predictor point is too far away from the curve, so the step width should be reduced (and vice versa)
                % However, this way of calculating r_pre has been removed since it showed similar or worse perfomance than step control method 'angle'
                % The remaining method is based on keeping norm(obj.p_y1-obj.yp) / norm(obj.dy0*obj.step_width) = norm_corr / obj.step_width constant, which leads to r_pre = const. / norm_corr
        
                % Step control parameters
                it_max = obj.step_control_param(1);                     % Maximal number of corrector iterations
                norm_corr_nom = obj.step_control_param(2);              % Nominal value of norm "corrector-vector"

                % Additional constraint: Number of corrector iterations must not exceed nominal number of corrector iterations
                if obj.p_it > it_max                                    % If number of corrector iterations exceeds nominal number of corrector iterations ...
                    r_pre = r_limit(1);                                 % Minimal factor of r is taken to reduce step width
                    display_message_info = 'too_much_iterations';       % Display message will be set below since message can appear for multiple step control methods
        
                % Main criterion: adapt step width based on keeping norm_corr / obj.step_width constant
                else
                    norm_corr = norm(obj.p_y1-obj.yp);                  % Norm of "corrector-vector"
                    r_pre = norm_corr_nom / norm_corr;                  % Calculate preliminary value of r by setting a nominal value for norm "corrector-vector". This leads to step widths being more regular
                    display_message_info = append('norm_corr = ', num2str(norm_corr), ', corr_it = ', num2str(obj.p_it));   % Set display message
                end
                d = 0;                                                  % Set damping factor (generally, no damping is required when using step control method 'norm_corrector')

           
            case 'angle'
                % Step control using the angle between the direction vectors at the last two points
                % If the angle is large, the discretization of the curve is too coarse with respect to the local curvature, so the step width should be reduced (and vice versa)
                
                % Step control parameters
                it_max = obj.step_control_param(1);                      % Maximal number of corrector iterations
                alpha_nom = obj.step_control_param(2);                   % Nominal angle

                % Additional constraint: Number of corrector iterations must not exceed nominal number of corrector iterations
                if obj.p_it > it_max                                    % If number of corrector iterations exceeds nominal number of corrector iterations ...
                    r_pre = r_limit(1);                                 % Minimal factor of r is taken to reduce step width
                    display_message_info = 'too_much_iterations';       % Display message will be set below since message can appear for multiple step control methods

                % Main criterion: adapt step width using angle between the direction vectors at the last two points on curve (y0 and y0_old)
                else
                    alpha = acos( obj.dy0' * obj.p_dy_old );            % Angle between direction vectors at y0 and y0_old
                    if alpha == 0                                       % Needs to be checked in order to avoid division by 0
                        r_pre = r_limit(2);                             % Preliminary value of r can be set to maximum value if alpha = 0 (i.e. direction vectors dy0 and dy_old are equal)
                    else
                        r_pre = alpha_nom / alpha;                      % Calculate preliminary value of r using asymptotic expansion of alpha. alpha_nom sets the nominal (desired) angle
                    end
                    display_message_info = append('alpha = ', num2str(alpha*180/pi), ', corr_it = ', num2str(obj.p_it));    % Set display message
                end
                d = 0;                                                  % Set damping factor (generally, no damping is required when using step control method 'angle')


            case 'combination'                  
                % Step control using the number of corrector iterations and the angle between the direction vectors at the last two points
        
                % Step control parameters
                it_nom = obj.step_control_param(1);                      % Nominal number of corrector iterations
                alpha_nom = obj.step_control_param(2);                   % Nominal angle

                r_it = it_nom / obj.p_it;                               % Ratio of nominal number of corrector iterations to actual number of corrector iterations       
                alpha = acos( obj.dy0' * obj.p_dy_old );                % Angle between direction vectors at y0 and y0_old
                if alpha == 0                                           % Needs to be checked in order to avoid division by 0
                    r_alpha = r_limit(2);                               % r_alpha can be set to maximum value if alpha = 0 (i.e. direction vectors dy0 and p_dy_old are equal)
                else
                    r_alpha = alpha_nom / alpha;                        % Calculate preliminary value of r using asymptotic expansion of alpha. alpha_nom sets the nominal (desired) angle
                end
                if (r_it >= 1) && (r_alpha >= 1)                        % If both of the calculated r-values are >= 1 ...
                    r_pre = max([r_it,r_alpha]);                        % Take the r-value which is larger to ensure faster increase of step width
                else                                                    % In all other cases ...
                    r_pre = min([r_it,r_alpha]);                        % Take the r-value which is smaller to be on the safe side
                end
                d = 0.25;                                               % Set damping factor (default: 0.25)
                display_message_info = append('r_it = ', num2str(r_it), ', r_alpha = ', num2str(r_alpha), ', corr_it = ', num2str(obj.p_it));       % Set display message
            

            case 'pid'
                % Step control using PID control to adapt step width
                % This method was taken from  Valli, Elias, Carey, Coutinho (2009): "PID adaptive control of incremental and arclength continuation" and has been adapted to fit CoSTAR
                % ATTENTION: This method does not work well currently. The optimal values of the four parameters are still unknown

                % Step control parameters
                it_max = obj.step_control_param(1);                      % Maximal number of corrector iterations
                tol = obj.step_control_param(2);                         % Reference value
                kP = obj.step_control_param(3);   kI = obj.step_control_param(4);   kD = obj.step_control_param(5);

                % Additional constraint: Number of corrector iterations must not exceed nominal number of corrector iterations
                if obj.p_it > it_max                                    % If number of corrector iterations exceeds nominal number of corrector iterations ...
                    r_pre = r_limit(1);                                 % Minimal factor of r is taken to reduce step width
                    display_message_info = 'too_much_iterations';       % Display message will be set below since message can appear for multiple step control methods
                
                % Main criterion: adapt step width using a PID control algorithm
                else
                    % Calculate parameter e using the direction vectors at current point (y0) and previous point
                    obj.p_e = norm(obj.dy0 - obj.p_dy_old) / tol;               % norm(obj.dy0 - obj.p_dy_old) / norm(obj.dy0) = norm(obj.dy0 - obj.p_dy_old)
                
                    % Calculate preliminary value of r (r_pre)
                    if obj.p_e ~= 0                                             % If parameter e is unequal zero, r_pre can be calculated
                        prop = (obj.p_e_old/obj.p_e)^kP;                        % "Present"
                        int = (1/obj.p_e)^kI;                                   % "Past"
                        diff = (obj.p_e_old^2/(obj.p_e*obj.p_e_old_old))^kD;    % "Future" --> swap e_old_old and e_old in order to make sense?
                        r_pre = prop * int * diff;
                        display_message_info = append('e = ', num2str(obj.p_e), ', prop = ', num2str(prop), ', int = ', num2str(int), ', diff = ', num2str(diff), '. corr_it = ', num2str(obj.p_it));   % Set display message
                    else
                        r_pre = r_limit(2);                                     % If parameter e is zero (direction vectors dy0 and dy_old are equal), r_pre can be set to maximum value
                        display_message_info = append('e = 0, corr_it = ', num2str(obj.p_it));      % Set display message
                    end
                end
                d = 0;                                                  % Set damping factor (generally, no damping is required when using step control method 'pid')

        end


        % Calculate final factor r (obj.p_r) to adapt step width
        r_pre_limit = max( min(r_pre,r_limit(2)), r_limit(1) );         % If r_pre exceeds upper or lower limit, corresponding limit is taken. Otherwise, r_pre is taken
        obj.p_r = (1-d) * r_pre_limit + d * obj.p_r_old;                % Calculate final factor r. d is a damping factor controlling how much of the preceding r-value (obj.p_r_old) is taken into account

        % Calculate new step width
        step_width_pre = obj.step_width * obj.p_r;                                                          % New preliminary step width
        obj.step_width = max( min(step_width_pre,obj.step_width_limit(2)), obj.step_width_limit(1) );       % Set step_width. If step_width_pre exceeds limits, corresponding limit is taken. Otherwise, step_width_pre is taken


        % Display information
        if strcmpi(obj.display,'step_control_info')                     % Display information if message is desired

            if strcmpi(display_message_info, 'too_much_iterations')     % Set info message if corrector needed more iterations than desired
                display_message_info = append('Corrector needed more iterations than desired: corr_it = ', num2str(obj.p_it), ' > ', num2str(it_max), ' = corr_it_max');
            end

            if step_width_old == obj.step_width                         % If step width has not been adapted ...
                if obj.step_width == obj.step_width_limit(1)            % ... because bottom step width limit has been reached
                    disp(append('Step width not adapted since bottom step width limit has been reached (r_pre = ', num2str(r_pre), '). ', display_message_info))
                elseif obj.step_width == obj.step_width_limit(2)        % ... because upper step width limit has been reached
                    disp(append('Step width not adapted since upper step width limit has been reached (r_pre = ', num2str(r_pre), '). ', display_message_info))
                else                                                    % ... because factor r = 1
                    disp(append('Step width not adapted since step width is optimal (r_pre = ', num2str(r_pre), '). ', display_message_info))
                end

            else                                                        % If step width has been adapted
                disp(append('Step width adapted to step_width = ', num2str(obj.step_width), ' (r = ', num2str(obj.p_r), ', r_pre = ', num2str(r_pre), '). ', display_message_info));

            end

        end


end

end