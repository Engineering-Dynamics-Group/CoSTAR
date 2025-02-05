% Function calculates initial curve point
%
% @DYN: DynamicalSystem object
% @S:   Solution object
% @AM:  ApproxMethod object
% @ST:  Stability object

function [S,AM,DYN] = initial_solution(DYN,S,AM,ST)

if ~strcmpi(DYN.display,'off')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%                                 CoSTAR                                  %')
    disp('%              Continuation of Solution Torus AppRoximations              %')
    disp(append('%                              Version ',DYN.costar_version,'                                %'))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' ')
end


% Set function and options for fsolve
newtonOpts = optimoptions('fsolve','Display','iter-detailed','MaxFunEvals',1e5,'FiniteDifferenceType','forward');
y0 = [AM.iv;DYN.param{DYN.act_param}];

if strcmpi(DYN.approx_method,'shooting')
    Fcn = @(y) AM.fun_Jac_wrapper_init(y,y0);                           % Function wrapper for initial solution, if Jacobian is supplied
    newtonOpts.SpecifyObjectiveGradient = true;                         % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended

elseif strcmpi(DYN.approx_method,'finite-difference')                   % Special corrector function due to specification of Jacobian matrix
    newtonOpts.SpecifyObjectiveGradient = true;                         % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended
    Fcn = @(y) AM.corr_fun_init_FDM(y,y0);

else
    Fcn = @(y)[AM.res(y);y(end)-y0(end)];                               % Last entry is not necessary, but this way the Jacobian has the correct dimension

end


%%%%%%%%%%%%%%%%%%%%%%%%%% Corrector %%%%%%%%%%%%%%%%%%%%%%%%%%

warn_msg = cell(0,0);                                                   % Initialize

write_log(DYN,'initialize')                                             % Initialize log file
if strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final')
    [fsolve_text,y,~,newton_flag,~,J] = evalc('fsolve(Fcn,y0,newtonOpts)');
    write_log(DYN,'fsolve is trying to find an initial solution ...')
    write_log(DYN,fsolve_text(1:end-1))
else
    write_log(DYN,'diary_on')                                           % Start recording of command window for log file
    disp('fsolve is trying to find an initial solution ...')            % This text is automatically printed to the log file due to diary function
    [y,~,newton_flag,~,J] = fsolve(Fcn,y0,newtonOpts);
    write_log(DYN,'diary_off')                                          % Stop recording of command window for log file
end
% checkGradients_opts = optimoptions('fsolve',FiniteDifferenceType='forward'); checkGradients(Fcn,y,checkGradients_opts,Display='on',Tolerance=1e-5);   % Check the user-defined Jacobian matrix

% No initial solution found or solution found but Jacobian can be undefined
if (newton_flag < 1) || (newton_flag == 2)
    if newton_flag < 1
        error_text = 'ERROR: No initial solution found!';                       % Set error text
        stopping_msg = 'CoSTAR stopped because corrector did not converge.';    % Set stopping message
    elseif newton_flag == 2
        error_text = 'ERROR: Equation solved, but change in y smaller than the specified tolerance, or Jacobian at y is undefined!';    % Set error text
        stopping_msg = 'CoSTAR stopped because Jacobian can be undefined at initial solution.';                                         % Set stopping message
    end
    write_log(DYN,error_text)                                           % Write error text in log file
    write_log(DYN,'finalize_error',stopping_msg)                        % Finalize log file with error message
    if ~strcmpi(DYN.display,'off')
        disp(error_text); disp(' '); disp(stopping_msg); disp(' ')      % Display error text and message
        disp('-----------------------------------------------------')
        disp('--------------  Finished with error!  ---------------');
        disp('-----------------------------------------------------')
        disp(' ');
    end
    error(append(error_text(8:end),' ',stopping_msg))                   % Throw error

% Initial solution found, but fsolve exitflag is not default exitflag
elseif newton_flag == 3
    warn_msg{end+1} = 'WARNING: Equation solved, but change in residual is smaller than specified tolerance!';
    S.warnings{end+1} = warn_msg{end}(10:end);                          % Save warning in Solution object
elseif newton_flag == 4
    warn_msg{end+1} = 'WARNING: Equation solved, but magnitude of search direction is smaller than specified tolerance!';
    S.warnings{end+1} = warn_msg{end}(10:end);                          % Save warning in Solution object
end

info_text = 'Initial solution found!';                                  % This message is updated if stability is successfully computed


%%%%%%%%%%%%%%%%%%%%%%% Frequency Check %%%%%%%%%%%%%%%%%%%%%%%

stopping_msg = check_freq(DYN,y);
if ~isempty(stopping_msg)                                               % If frequency(s) are smaller than frequency limit
    DYN.cont = 'off';                                                   % Deactivate continuation (if it was set to 'on')
    warn_msg{end+1} = 'WARNING: Small or negative frequency(s) detected for initial solution!';      % Set warning message
    S.warnings{end+1} = warn_msg{end}(10:end);                          % Save warning in Solution object
    % Everything else (displaying warning, writing log file etc.) has to be done after this if...else... block

else


    %%%%%%%%%%%%%%%%%%%%%%%%%% Error Control Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
    % This loop is slightly different from the loop in the Continuation class
    % since it is does not update the dimension of the last curve point

    added_output = cell(0,0);                                               % Initialize (needed for saving)

    if strcmpi(AM.error_control,'on')

        iterate = 1;
        counter = 0;
        increase = 1;           % These values help to avoid a constant switching between de- and increasing, if error tolerance are set too close to each other
        decrease = 1;           % If we start increasing or decreasing once, the other one is not possible anymore

        while iterate

            counter = counter + 1;
            err = AM.IF_estimate_error(y,DYN);
            err_control_text = append('Current Error: ',num2str(err),' -- Number of Harmonics: ',num2str(size(AM.hmatrix,2)-1));
            if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                disp(err_control_text);
            end
            write_log(DYN,err_control_text)

            if (err > AM.error_limit(2)) && (counter < AM.ec_iter_max+1) && increase
                err_control_text = append('Increase Discretization (Iteration: ',num2str(counter),')');
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(err_control_text);
                end
                write_log(DYN,err_control_text)
                [y0,iterate] =  AM.IF_increase_discretization(y,DYN); %this might also update properties in AM
                decrease = 0;
            elseif (err < AM.error_limit(1)) && (counter < AM.ec_iter_max+1) && decrease
                err_control_text = append('Decrease Discretization (Iteration: ',num2str(counter),')');
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(err_control_text);
                end
                write_log(DYN,err_control_text)
                [y0,iterate] = AM.IF_decrease_discretization(y,DYN); %this might also update properties in AM
                increase = 0;
            else
                iterate = 0;
            end

            if iterate          % Recalculate a new curve point with de-/increased discretisation
                % If autonomous y_old contains predicted solution
                if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
                    AM.IF_up_res_data(obj,DYN);                     % !!! Check this line if error_control becomes available for quasiperiodic shooting !!!
                    AM.y_old = y0;
                    Fcn = @(y)AM.fun_Jac_wrapper_init(y,y0,DYN);    % Function wrapper for initial solution, if Jacobian is supplied
                    newtonOpts.SpecifyObjectiveGradient = true;
                else
                    AM.iv = y0(1:(end-1));                          % Archive initial solution
                    Fcn = @(y)[AM.res(y);y(end)-y0(end)];           % Define corrector-function containing the residual function and the subspace-constraint
                end

                newtonOpts.Display = 'off';                         % Deactivate fsolve output
                [y_er,~,newton_flag_er,~,J_er] = fsolve(Fcn,y0,newtonOpts);  % Solve corrector-function

                % This if...else below does not work yet!
                % Reason: AM.iv makes changes to AM object. When fsolve fails and the new solution is not accepted, we need to revert the mentioned changes.
                %         If we do not revert the changes, there will be an error in stability computation! But how do we revert the changes?
                % if (newton_flag_er > 0) && (newton_flag_er ~= 2)  % Accept the result from fsolve if exitflag is 1, 3 or 4
                    y = y_er;
                    newton_flag = newton_flag_er;
                    J = J_er;
                % else                                                % If exitflag from fsolve is < 0 or = 2
                %     warn_msg{end+1} = 'WARNING: Error control stopped early or failed for initial solution!';
                %     S.warnings{end+1} = warn_msg{end}(10:end);      % Save warning in Solution object
                %     break                                           % Immediately break the loop. The last accepted solution is returned
                % end

            end

        end

        if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
            disp(' ');
        end
        write_log(DYN,'')
        added_output{1} = err;

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%% STABILITY %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Has to be called after the error_control step

    if strcmpi(DYN.stability,'on')

        [multipliers,vectors,n_unstable,stability_flag] = ST.calc_stability(y,J);
        added_output{2} = multipliers;
        added_output{3} = n_unstable;
        added_output{4} = stability_flag;
        added_output{5} = vectors;

        if stability_flag == 0                                                              % Stability computation failed
            warn_msg{end+1} = 'WARNING: Stability computation failed for the initial solution!';    % Set warning message
            S.warnings{end+1} = warn_msg{end}(10:end);                                      % Save warning in Solution object
        elseif n_unstable == 0                                                              % Stable solution found
            info_text = 'Initial solution (stable) found!';                                 % Update info text
        elseif n_unstable > 0                                                               % Unstable solution found
            info_text = 'Initial solution (unstable) found!';                               % Update info text
        end

    end


    %%%%%%%%%%%%%%%%%%%% Save initial solution %%%%%%%%%%%%%%%%%%%%

    if strcmpi(DYN.sol_type,'quasiperiodic') && strcmpi(DYN.approx_method,'shooting')
        S.IF_arch_init_data(y,J,AM.Ik,newton_flag,AM.phi,DYN,added_output);
    else
        S.IF_arch_init_data(y,J,newton_flag,DYN,AM,added_output);
    end

end


%%%%%%%%%%% Display info and warning and write log %%%%%%%%%%%%

if ~strcmpi(DYN.display,'off')
    disp(info_text); disp(' ');                 % Display info text
end
write_log(DYN,info_text)                        % Write info text to log
if ~isempty(warn_msg)
    for i = 1:numel(warn_msg)
        warning(warn_msg{i}(10:end)); disp(' ') % Display warning
        write_log(DYN,append('\n',warn_msg{i})) % Write warning in log file
    end
end


%%%%%%%%%%%%%%%%% Stopping (no continuation) %%%%%%%%%%%%%%%%%%

if strcmpi(DYN.cont,'off')
    if isempty(stopping_msg)                    % Set stopping message if it has not already been set above
        stopping_msg = 'CoSTAR stopped after initial solution because continuation was turned off.';  
    end
    write_log(DYN,'finalize',stopping_msg)      % Finalize log file with stopping message
    S.stopping_flag = stopping_msg;             % Save stopping message
    if ~strcmpi(DYN.display,'off')
        disp(stopping_msg)                      % Display stopping message
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

end