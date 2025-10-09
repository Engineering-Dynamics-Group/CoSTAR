% This function controls the error of the approximation by adapting the discretisation and reiterating the solution
%
% @obj:  Continuation class object
% @S:    Solution class object
% @AM:   ApproxMethod class object
% @DYN:  Dynamical System class object

function obj = error_control(obj,S,AM,DYN)

    iterate = 1;
    counter = 0;
    increase = 1;       % These values help to avoid a constant switching between de- and increasing, if error tolerance are set too close to each other
    decrease = 1;       % If we start increasing or decreasing once, the other one is not possible anymore

    while iterate

        counter = counter + 1;
        obj.IF_ec_save_reset_props(AM,'save');                  % Save the properties that are modified by the error control
        obj.p_error = AM.IF_estimate_error(obj.p_y1,DYN);
        err_control_text = append('Current Error: ',num2str(obj.p_error),' -- Number of Harmonics: ',num2str(size(AM.hmatrix,2)-1));
        write_log(DYN,err_control_text)
        if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
            disp(err_control_text); 
        end

        %%%%%%%%%% Error is too high %%%%%%%%%%
        if (obj.p_error > AM.error_limit(2)) && (counter < AM.ec_iter_max+1) && increase
            [obj.yp,iterate] = AM.IF_increase_discretization(obj.p_y1,DYN);                 % This might also update properties in AM
            [obj.y0,obj.dy0] = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.y0,obj.dy0);      % Updates the older data points to match the dimension of the new solution vector (depends on solution type)
            decrease = 0;
            err_control_text = append('Increase Discretization (Iteration: ',num2str(counter),')');
            write_log(DYN,err_control_text)
            if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                disp(err_control_text); 
            end

        %%%%%%%%%% Error is too low %%%%%%%%%%
        elseif (obj.p_error < AM.error_limit(1)) && (counter < AM.ec_iter_max+1) && decrease
            [obj.yp,iterate] = AM.IF_decrease_discretization(obj.p_y1,DYN);                 % This might also update properties in AM 
            [obj.y0,obj.dy0] = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.y0,obj.dy0);      % Updates the older data points to match the dimension of the new solution vector (depends on solution type)
            increase = 0;
            err_control_text = append('Decrease Discretization (Iteration: ',num2str(counter),')');
            write_log(DYN,err_control_text)
            if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                disp(err_control_text); 
            end 

        %%%%%%%%%% Error is fine %%%%%%%%%%
        else
            
            obj.p_ec_flag = 1;  % Set exit flag of error control
            break               % Stop the loop of adapting the discretisation (iterate = 0) and skip the following recalculation

        end

        
        %%%%%%%%%% Recompute solution with in-/decreased discretization %%%%%%%%%%
        if iterate == 1
            % Recalculate a new curve point with de-/increased discretisation if error is not fine
            if strcmpi(DYN.approx_method,'shooting')                        % Special corrector function for SHM due to specification of Jacobian matrix
                AM.IF_up_res_data(obj,DYN);                                 % Pass information to the ApproxMethod object
                Fcn = @(y) AM.fun_Jac_wrapper(y,obj);                       % Set functionwrapper to provide Jacobian
            elseif strcmpi(DYN.approx_method,'finite-difference')           % Special corrector function for FDM due to specification of Jacobian matrix
                AM.IF_up_res_data(obj);                                     % Pass information to the ApproxMethod object
                Fcn = @(y) AM.corr_fun_FDM(y,obj);                          % Set corrector-function
            else
                AM.IF_up_res_data(obj);                                     % Pass information to the ApproxMethod object
                Fcn = @(y)[AM.res(y);obj.sub_con(y,obj)];                   % Define corrector-function containing the residual function and the subspace-constraint
            end

            [y1,~,exit_flag,output,J1] = fsolve(Fcn,obj.yp,obj.fsolve_opts);          % Solve corrector-function

            % Exit flag handling I
            if any(exit_flag == [1,3,4])                                    % Accept the result from fsolve if exitflag is 1, 3 or 4
                obj.p_y1 = y1;
                obj.p_newton_flag = exit_flag;
                obj.p_output = output;
                obj.p_J1 = J1;
                obj.p_y0_old{2} = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.p_y0_old{2});  % We must update p_y0_old to the new solution dimension
                obj.p_y0_old{3} = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.p_y0_old{3});  % We must update p_y0_old to the new solution dimension
                obj.p_ec_flag = 1;
                if exit_flag == 3
                    warn_text = append('Equation solved for new discretization in error control, but change in residual is smaller than specified tolerance (fsolve exit_flag = 3)!');
                    write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
                    S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
                    obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
                    warning(warn_text);                                                                         % Display warning
                elseif exit_flag == 4
                    warn_text = append('Equation solved for new discretization in error control, but magnitude of search direction is smaller than specified tolerance (fsolve exit_flag = 4)!');
                    write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
                    S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
                    obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
                    warning(warn_text);                                                                         % Display warning
                end

            % Exit flag handling II: If exitflag from fsolve is <= 0 or = 2: Issue warning, reduce step width (if possible) and terminate error control
            else
                warn_text = append('Error control failed for solution Iter = ',num2str(obj.p_local_cont_counter+1),', because fsolve returned exit_flag = ',num2str(exit_flag),'!');
                write_log(DYN,append('WARNING: ',warn_text))                % Write warning in log file
                S.warnings{end+1} = warn_text;                              % Save warning in Solution object
                obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
                warning(warn_text)                                          % Display warning
            end
               

        %%%%%%%%%% In-/Decrease of discretization failed %%%%%%%%%%
        elseif iterate == 0
            warn_text = append('In-/Decrease of discretization failed during error control for solution Iter = ',num2str(obj.p_local_cont_counter+1),'!');
            write_log(DYN,append('WARNING: ',warn_text))                    % Write warning in log file
            S.warnings{end+1} = warn_text;                                  % Save warning in Solution object
            obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
            warning(warn_text)                                              % Display warning
            
            
        %%%%%%%%%% Decreasing not possible: Minimum discretization reached %%%%%%%%%%
        elseif iterate == -1
            info_text = 'Number of harmonics cannot be reduced any further since minimum number has been reached.';        % Set info text
            write_log(DYN,info_text)                                        % Write info text in log file
            if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                disp(info_text);                                            % Display information
            end
            obj.p_ec_flag = 1;  % Set exit flag of error control
            break               % We can return to m_continuation now and can skip Exit flag handling III

        end


        %%%%%%%%%% Exit flag handling III: Reduce step width or terminate %%%%%%%%%%
        if iterate == 0 || ~any(exit_flag == [1,3,4])

            % Reduce step width if possible - else: stop continuation
            if obj.step_width > obj.step_width_limit(1,1)
                obj.IF_ec_save_reset_props(AM,'reset');                                                     % Reset the modified properties -> needed to continue since new solution is not accepted
                step_width_pre = 0.5.*obj.step_width;                                                       % New preliminary step width
                obj.step_width = max([step_width_pre,obj.step_width_limit(1)]);                             % Set step_width. If new preliminary step width falls below minimal step width, take minimal step width
                obj.p_convergence = 0;                                                                      % Set property p_convergence to zero (for resetting the step_width after convergence)
                info_text = append('Step width adapted to step_width = ',num2str(obj.step_width),', because error control failed! Computing a new solution with reduced step width.');
                write_log(DYN,info_text)                                                                    %write info text in log file
                if strcmpi(DYN.display,'step-control') || strcmpi(DYN.display,'full'); disp(info_text); end %display info text
                obj.p_ec_flag = 0;                                                                          % Set exit flag of error control
            else
                stopping_msg = 'CoSTAR stopped because error control failed and step width has reached minimal value.';   %set stopping message
                obj.p_stopping_flag = stopping_msg;                                                         %save stopping message
                if ~strcmpi(DYN.display,'off'); disp(' '); disp(stopping_msg); end                          %display stopping message
                write_log(DYN,'finalize',stopping_msg)                                                      %finalize log file with warning message and message
                obj.p_ec_flag = -1;                                                                         % Set exit flag of error control
            end

        end

    end

end