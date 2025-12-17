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
if(~strcmpi(obj.pred,'tangent'))
    obj = obj.initial_slope(DYN,AM);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  CONTINUATION  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while  obj.p_contDo
        
    %%%%%%%%%%%%%  PREDICTOR AND STEP CONTROL  %%%%%%%%%%%%
    if any(obj.p_newton_flag == [0,1,3,4]) && obj.p_ec_flag == 1    %direction vector needs to be calculated in the first loop and every time a new solution has been found
        obj.direction_vector();                 %calculate direction vector
    end
   
    if any(obj.p_newton_flag == [1,3,4]) && obj.p_ec_flag == 1     %stepcontrol may be called if corrector converged and error control is fine (initialised: obj.p_newton_flag = 0, obj.p_ec_flag = 1)
        obj.stepcontrol(DYN);                   %adapt step width
        obj.p_convergence = 1;                  %reset convergence property if corrector did not converge previously
    end     
    
    obj.predictor();                            %calculate predicted point

    stopping_msg = check_freq(DYN,obj.yp);                                      %check the frequencies at the predictor point (not done in predictor to be able to break the while loop)
    if ~isempty(stopping_msg)                                                   %if frequency(s) are smaller than frequency limit
        warn_msg = append('WARNING: Small or negative frequency(s) detected for predictor point after Iter = ',num2str(obj.p_local_cont_counter),'!');
        S.warnings{end+1} = warn_msg(10:end);                                   %save warning in Solution object
        obj.p_stopping_flag = stopping_msg;                                     %save stopping message
        disp(' '); warning(warn_msg(10:end));                                   %display warning
        if ~strcmpi(DYN.display,'off'); disp(' '); disp(stopping_msg); end      %display stopping message
        write_log(DYN,'finalize',append(warn_msg,'\n\n',stopping_msg))          %finalize log file with warning message and message
        break                                                                   %break the continuation loop (any further computation below might fail)
    end


    %%%%%%%%%%%%%%%%%%%%%%  CORRECTOR  %%%%%%%%%%%%%%%%%%%%
    if strcmpi(DYN.approx_method,'shooting')                                    %special corrector function for quasi-periodic shooting
        obj.fsolve_opts.MaxIter = 50;
        AM.IF_up_res_data(obj,DYN);                                             %pass information to the ApproxMethod object
        Fcn = @(y) AM.fun_Jac_wrapper(y,obj);                                   %set functionwrapper to provide Jacobian
        obj.fsolve_opts.SpecifyObjectiveGradient = true;                        %Jacobian matrix is passed by the user
    elseif strcmpi(DYN.approx_method,'finite-difference')                       %special corrector function for FDM due to specification of Jacobian matrix
        AM.IF_up_res_data(obj);                                                 %pass information to the ApproxMethod object
        obj.fsolve_opts.SpecifyObjectiveGradient = true;                        %Jacobian matrix is passed by the user
        Fcn = @(y) AM.corr_fun_FDM(y,obj);                                      %set corrector-function
    else
        AM.IF_up_res_data(obj);                                                 %pass information to the ApproxMethod object
        Fcn = @(y)[AM.res(y);obj.sub_con(y,obj)];                               %define corrector-function containing the residual function and the subspace-constraint
    end
    
    [obj.p_y1,~,obj.p_newton_flag,obj.p_output,obj.p_J1] = fsolve(Fcn,obj.yp,obj.fsolve_opts);      %solve corrector function

    
    %%%%%%%%%%%%  EXITFLAG < 1 OR EXITFLAG = 2  %%%%%%%%%%%
    if ((obj.p_newton_flag < 1) || (obj.p_newton_flag == 2)) && (obj.step_width <= obj.step_width_limit(1,1))               %if step width is already <= minimal step width 
        if obj.p_newton_flag < 1
            warn_msg = append('WARNING: No solution found for Iter = ',num2str(obj.p_local_cont_counter+1),' (fsolve exit_flag = ',num2str(obj.p_newton_flag),')!');
            stopping_msg = 'CoSTAR stopped because corrector did not converge and step width has reached minimal value.';   %set stopping message
        elseif obj.p_newton_flag == 2
            warn_msg = append(['WARNING: Equation solved for Iter = ',num2str(obj.p_local_cont_counter+1),', but ' ...      %set warning message
                              'change in y smaller than the specified tolerance, or Jacobian at y is undefined (fsolve exit_flag = 2)!']);
            stopping_msg = 'CoSTAR stopped because Jacobian can be undefined and step width has reached minimal value.';    %set stopping message
        end
        S.warnings{end+1} = warn_msg(10:end);                                                       %save warning in Solution object
        obj.p_stopping_flag = stopping_msg;                                                         %save stopping message
        disp(' '); warning(warn_msg(10:end));                                                       %display warning
        if ~strcmpi(DYN.display,'off'); disp(' '); disp(stopping_msg); end                          %display stopping message
        write_log(DYN,'finalize',append(warn_msg,'\n\n',stopping_msg))                              %finalize log file with warning message and message
        break                                                                                       %immediately break while loop, because everything else can lead to errors
    
    elseif ((obj.p_newton_flag < 1) || (obj.p_newton_flag == 2)) && (obj.step_width > obj.step_width_limit(1,1))            %if fsolve did not converge and step width is above minimal step width 
        if obj.p_newton_flag < 1
            warn_text = append('No solution found for Iter = ',num2str(obj.p_local_cont_counter+1),' (fsolve exit_flag = ',num2str(obj.p_newton_flag),')!');
        elseif obj.p_newton_flag == 2
            warn_text = append(['Equation solved for Iter = ',num2str(obj.p_local_cont_counter+1),', but ' ...              %set warning message
                                'change in y smaller than the specified tolerance, or Jacobian at y is undefined (fsolve exit_flag = 2)!']);
        end
        write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
        S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
        obj.p_last_msg = sprintf('%s%s%s\n',obj.p_last_msg,append('Warning: ',warn_text),' ');      % Save the warning in the "last messages" property
        warning(warn_text);
        step_width_pre = 0.5.*obj.step_width;                                                       %new preliminary step width
        obj.step_width = max([step_width_pre,obj.step_width_limit(1)]);                             %set step_width. If new preliminary step width falls below minimal step width, take minimal step width
        obj.p_convergence = 0;                                                                      %set property p_convergence to zero (for resetting the step_width after convergence)
        info_text = append('Step width adapted to stepwidth = ',num2str(obj.step_width),', because corrector did not converge or Jacobian can be undefined!');
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
    stopping_msg = check_freq(DYN,obj.p_y1);                                    %check the frequencies of the solution
    if ~isempty(stopping_msg)                                                   %if frequency(s) are smaller than frequency limit
        warn_msg = append('WARNING: Small or negative frequency(s) detected for solution Iter = ',num2str(obj.p_local_cont_counter+1),'!');
        S.warnings{end+1} = warn_msg(10:end);                                   %save warning in Solution object
        obj.p_stopping_flag = stopping_msg;                                     %save stopping message
        disp(' '); warning(warn_msg(10:end));                                   %display warning
        if ~strcmpi(DYN.display,'off'); disp(' '); disp(stopping_msg); end      %display stopping message
        write_log(DYN,'finalize',append(warn_msg,'\n\n',stopping_msg))          %finalize log file with warning message and stopping message
        break                                                                   %break (stop) the continuation loop (any further computation below might fail)
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%  ERROR CONTROL  %%%%%%%%%%%%%%%%%%%%%%%
    %Control the error by adapting the discretisation. If the ansatz function is adapted, error_control solves the equation system again
    %and computes new values, which get saved and iterated, which is why the error_control must be called before IF_arch_data and iterate_data.
    if strcmpi(AM.error_control,'on')
        obj.error_control(S,AM,DYN);
        if obj.p_ec_flag == 0                                                   %Error control failed and step size is larger than minimum: Reduce step size and try again (start a new iteration)
            continue                                                            %Skip the following methods and start the next loop (try again with reduced step width)
        elseif obj.p_ec_flag == -1                                              %Error control failed and step size is minimal: Stop continuation
            break                                                               %Skip the following methods because the error of the solution is too high (i.e. do not accept the solution)
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%  STABILITY  %%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate the Lyapunov stability. If a bifurcation point is found, the location of the point can be iterated
    %IMPORTANT: Has to be called after the error_control step
    if strcmpi(DYN.stability,'on'); obj.bifurcation_stability(DYN,AM,S,ST); end


    %%%%%%%%%%%%%%%%%%%%%%%%  SYNCHRONIZATION  %%%%%%%%%%%%%%%%%%%%%%
    %If quasi-periodic solution is continued, check if solution
    %synchronizes
    if strcmpi(DYN.synchronization,'on'); obj.checkSync(DYN,AM,S); end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  SAVING  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    S.IF_arch_data(obj,DYN,AM);                                                 %store calculated data point in Solution object


    %%%%%%%%%%  PLOTTING, CHECKING LIMITS AND DISPLAY INFO  %%%%%%%%%
    if strcmpi(obj.plot,'on'); obj.plot_contplot(S,DYN); end                    %display a continuation plot if desired

    obj.iterate_data();                                                         %store calculated data point as current data point for next continuation step

    obj.check_limits(DYN);                                                      %check limits for bifurcation parameter and number of steps and display information
    
end


%%%%%%%%%%%%%%%  END CONTINUATION  %%%%%%%%%%%%%
S.stopping_flag = obj.p_stopping_flag;          %save stopping message in Solution object

if ~strcmpi(DYN.display,'off')
    disp(' ')
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