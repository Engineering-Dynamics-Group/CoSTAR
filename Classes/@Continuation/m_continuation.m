% Main method which continues a curve defined by the roots of an algebraic
% system of equations. This function needs an object of Solution (S) and an
% object of ApproxMethod (AM). m_continuation returns an object of Solution
% which contains the continued curve
%
%@obj:      Continuation class object
%@DYN:      DynamicalSystem object
%@S:        Solution subclass object
%@AM:       ApproxMethod subclass object
%@ST:       Stability subclass object
%
%@S:        Solution subclass object

function  S = m_continuation(obj,DYN,S,AM,ST)


%%%%%%%%%%%%%%%%  INITIALISATION  %%%%%%%%%%%%%%
obj.y0  = S.y0;                                                         %Get active curve point from object of Solution class
obj.p_y0_old = num2cell([zeros(numel(obj.y0),2),obj.y0],1);             %Must be initialized as 1x3 cell array for the error control to work
obj.p_mu0 = S.y0(end,1);
if isa(S.J,'cell'); obj.p_J0  = S.J{1,1}; else; obj.p_J0  = S.J; end
if strcmpi(DYN.stability,'on') 
    obj.p_stability_flag_old = S.stability_flag;
    if strcmpi(ST.iterate_bfp,'on')
        obj.p_n_unstable_0 = S.n_unstable; 
        ST.update_curve_container(DYN,AM,S.arclength,obj.y0,obj.y0,S.multipliers,obj.p_n_unstable_0); %Update the curve container with the first point
    end
end  

cont_header = sprintf('\n%s\n%s\n%s\n',...
                      '---------------------------------------------------------------------------',...
                      '-----------------------------  Continuation  ------------------------------',...
                      '---------------------------------------------------------------------------');
write_log(DYN,cont_header)
if ~strcmpi(DYN.display,'off')
    disp('-----------------------------------------------------')
    disp('------------------  Continuation  -------------------')
    disp('-----------------------------------------------------')
    disp(' ')
end

% Calculate second curve point with differential perturbation dmu to
% calculate initial slope for secant predictor
if ~strcmpi(obj.predictor,'tangent')
    obj = obj.initial_slope(DYN,AM);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  CONTINUATION  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while  obj.p_contDo
        
    %%%%%%%%%%%%%  PREDICTOR AND STEP CONTROL  %%%%%%%%%%%%
    if any(obj.p_newton_flag == [1,3,4]) && (obj.p_error_flag == 1)    %direction vector needs to be calculated in the first loop and every time a new solution has been found
        obj.direction_vector();                 %calculate direction vector
    end
   
    if any(obj.p_newton_flag == [1,3,4]) && (obj.p_error_flag == 1) && (obj.p_local_cont_counter > 1)  %stepcontrol may be called if corrector converged and error control is fine, but not in the first loop
        obj.stepcontrol(DYN);                   %adapt step width
        obj.p_convergence = 1;                  %reset convergence property if corrector did not converge previously
    end     
    
    obj.predictor_point();                      %calculate predictor point

    [obj.p_stopping_flag,obj.p_error_flag] = check_freq(DYN,obj.yp);            %check the frequencies at the predictor point (not done in predictor to be able to break the while loop)
    if obj.p_error_flag == 0                                                    %if frequency(s) are smaller than frequency limit
        obj.p_error_msg = append('ERROR: Small or negative frequency(s) detected for predictor point after Iter = ',num2str(obj.p_local_cont_counter),'!');
        break                                                                   %break the continuation loop (any further computation below might fail)
    end


    %%%%%%%%%%%%%%%%%%%%%%  CORRECTOR  %%%%%%%%%%%%%%%%%%%%
    AM.IF_up_res_data(obj,DYN);                                                 %update AM properties
    Fcn = @(y) AM.res_fun(y,obj);                                               %function wrapper to set the complete residuum function
    
    [obj.p_y1,~,obj.p_newton_flag,obj.p_output,obj.p_J1] = fsolve(Fcn,obj.yp,obj.fsolve_opts);      %solve corrector function

    
    %%%%%%%%%%%%  EXITFLAG < 1 OR EXITFLAG = 2  %%%%%%%%%%%
    if ((obj.p_newton_flag < 1) || (obj.p_newton_flag == 2)) && (obj.step_width <= obj.step_width_limit(1,1))               %if step width is already <= minimal step width 
        obj.p_error_flag = 0;                                                                                               %set error flag to indicate a critical error
        if obj.p_newton_flag < 1
            obj.p_error_msg = append('ERROR: No solution found for Iter = ',num2str(obj.p_local_cont_counter+1),' (fsolve exit_flag = ',num2str(obj.p_newton_flag),')!');
            obj.p_stopping_flag = 'CoSTAR stopped because corrector did not converge and step width has reached minimal value.';
        elseif obj.p_newton_flag == 2
            obj.p_error_msg = append(['ERROR: Equation solved for Iter = ',num2str(obj.p_local_cont_counter+1),', but ' ...
                                      'change in y smaller than the specified tolerance, or Jacobian at y is undefined (fsolve exit_flag = 2)!']);
            obj.p_stopping_flag = 'CoSTAR stopped because Jacobian can be undefined and step width has reached minimal value.';
        end
    
    elseif ((obj.p_newton_flag < 1) || (obj.p_newton_flag == 2)) && (obj.step_width > obj.step_width_limit(1,1))            %if fsolve did not converge and step width is above minimal step width 
        if obj.p_newton_flag < 1
            warn_text = append('No solution found for Iter = ',num2str(obj.p_local_cont_counter+1),' (fsolve exit_flag = ',num2str(obj.p_newton_flag),')! Trying again with reduced step width.');        
        elseif obj.p_newton_flag == 2
            warn_text = append(['Equation solved for Iter = ',num2str(obj.p_local_cont_counter+1),', but change in y is smaller than the specified ' ...
                                'tolerance, or Jacobian at y is undefined (fsolve exit_flag = 2)! Trying again with reduced step width.']);
        end
        write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
        S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
        obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
        warning(warn_text);
        step_width_pre = 0.5.*obj.step_width;                                                       %new preliminary step width
        obj.step_width = max([step_width_pre,obj.step_width_limit(1)]);                             %set step_width. If new preliminary step width falls below minimal step width, take minimal step width
        obj.p_convergence = 0;                                                                      %set property p_convergence to zero (for resetting the step_width after convergence)
        info_text = append('Step width adapted to stepwidth = ',num2str(obj.step_width),'.');
        write_log(DYN,info_text)                                                                    %write info text in log file
        if strcmpi(DYN.display,'step-control') || strcmpi(DYN.display,'full'); disp(info_text); end %display info text
        continue                                                                                    %skip the remaining code and start the next loop (try again with reduced step width)
    
    end


    %%%%%%%%%%%%%%%%%%  FSOLVE CONVERGED  %%%%%%%%%%%%%%%%%
    % A user-defined Jacobian matrix can be checked here by executing the next line (since R2023b: function checkGradients is recommended instead of obj.fsolve_opts.CheckGradients = true)
    % checkGradients_opts = optimoptions('fsolve',FiniteDifferenceType='forward'); checkGradients(Fcn,obj.p_y1,checkGradients_opts,Display='on',Tolerance=1e-5);
    if obj.p_newton_flag == 3
        warn_text = append('Equation solved for Iter = ',num2str(obj.p_local_cont_counter+1),', but change in residual is smaller than specified tolerance (fsolve exit_flag = 3)!');
        write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
        S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
        obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
        warning(warn_text);                                                                         % Display warning
    elseif obj.p_newton_flag == 4
        warn_text = append('Equation solved for Iter = ',num2str(obj.p_local_cont_counter+1),', but magnitude of search direction is smaller than specified tolerance (fsolve exit_flag = 4)!');
        write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
        S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
        obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
        warning(warn_text);                                                                         % Display warning
    end


    %IMPORTANT: order of calling the following methods must not be changed

    %%%%%%%%%%%%%%%%%%%%%%%  FREQUENCY CHECK  %%%%%%%%%%%%%%%%%%%%%%%
    if obj.p_error_flag == 1
        [obj.p_stopping_flag,obj.p_error_flag] = check_freq(DYN,obj.p_y1);              % Check the frequencies of the solution
        if obj.p_error_flag == 0                                                        % If frequency(s) are smaller than frequency limit
            obj.p_error_msg = append('ERROR: Small or negative frequency(s) detected for solution Iter = ',num2str(obj.p_local_cont_counter+1),'!');
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%  ERROR CONTROL  %%%%%%%%%%%%%%%%%%%%%%%
    %Control the error by adapting the discretisation. If the ansatz function is adapted, error_control solves the equation system again
    %and computes new values, which get saved and iterated, which is why the error_control must be called before IF_arch_data and iterate_data.
    if strcmpi(AM.error_control,'on'); obj.error_control(S,AM,DYN);
        if obj.p_error_flag == 2                                                %Error control failed and step size is larger than minimum: Reduce step size and try again (start a new iteration)
            continue                                                            %Skip the following methods and start the next loop (try again with reduced step width)
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%  STABILITY  %%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the Lyapunov stability. If a bifurcation point is found, the location of the point can be iterated
    %IMPORTANT: Has to be called after the error_control step
    if strcmpi(DYN.stability,'on')
        obj.bifurcation_stability(DYN,AM,S,ST);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK RETURN %%%%%%%%%%%%%%%%%%%%%%%%
    % Checks if curve returns to a previously computed part of the curve
    % IMPORTANT: Has to be called after error control and stability, but before saving
    if strcmpi(AM.error_control,'off') && (obj.p_error_flag == 1)                              
        obj.check_return(DYN,S);                                                %check if curve returned on itself
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SAVING  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    S.IF_arch_data(obj,DYN,AM);                                                 %store calculated data point in Solution object
    if strcmpi(DYN.save,'on')
        save(append('CoSTAR_Save_',DYN.DYN_id,'.mat'),'DYN','S');               %save DYN and S in a mat-file (overrides the existing file)
    end


    %%%%%%%%%%  PLOTTING, CHECKING LIMITS AND DISPLAY INFO  %%%%%%%%%
    if strcmpi(obj.plot,'on') && (obj.p_error_flag == 1)
        obj.plot_contplot(S,DYN);                                               %display a continuation plot if desired
    end

    obj.iterate_data();                                                         %store calculated data point as current data point for next continuation step

    obj.check_limits(DYN);                                                      %check limits for bifurcation parameter and number of steps and display information


end


%%%%%%%%%%%%%%%  END CONTINUATION  %%%%%%%%%%%%%
S.stopping_flag = obj.p_stopping_flag;              % Save stopping message in Solution object

if obj.p_error_flag == 0                            % Stopping in case of a critical error
    S.warnings{end+1} = obj.p_error_msg;            % Save error message in Solution object
    disp(' '); warning(obj.p_error_msg);            % Display error
    write_log(DYN,'finalize_error',append(obj.p_error_msg,'\n\n',obj.p_stopping_flag))  % Finalize log file with error message and stopping message
    if ~strcmpi(DYN.display,'off')
        disp(' '); disp(obj.p_stopping_flag); disp(' ')    % Display error text and stopping message
        disp('-----------------------------------------------------')
        disp('--------------  Finished with error!  ---------------')
        disp('-----------------------------------------------------')
        disp(' '); disp(' ')
    end

elseif obj.p_error_flag == 1
    write_log(DYN,'finalize',obj.p_stopping_flag)
    if ~strcmpi(DYN.display,'off')                  % Regular stopping
        disp(' '); disp(obj.p_stopping_flag); disp(' ')
        disp('-----------------------------------------------------')
        if ~isempty(lastwarn)
            disp('------------  Finished with warning(s)!  ------------');
        else
            disp('-------------- Successfully finished! ---------------')
        end
        disp('-----------------------------------------------------')
        disp(' '); disp(' ')
    end
end

if strcmpi(DYN.save,'on')
    save(append('CoSTAR_Save_',DYN.DYN_id,'.mat'),'DYN','S');               %save again to also catch the stopping flag saved in S
end


end