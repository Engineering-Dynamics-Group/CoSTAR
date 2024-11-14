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
    disp(append('%                             Version ',DYN.costar_version,'                               %'))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' ')
end


% Set function and options for fsolve
newtonOpts = optimoptions('fsolve','Display','iter-detailed','MaxFunEvals',1e5,'FiniteDifferenceType','forward');
y0 = [AM.iv;DYN.param{DYN.act_param}];

if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
    Fcn = @(y)AM.fun_Jac_wrapper_init(y,y0,DYN);                        % Function wrapper for initial solution, if Jacobian is supplied
    newtonOpts.SpecifyObjectiveGradient = true;

elseif strcmpi(DYN.approx_method,'finite-difference')                   % Special corrector function due to specification of Jacobian matrix
    newtonOpts.SpecifyObjectiveGradient = true;                         % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended
    Fcn = @(y) AM.corr_fun_init_FDM(y,y0);

else
    Fcn = @(y)[AM.res(y);y(end)-y0(end)];                               % Last entry is not necessary, but this way the Jacobian has the correct dimension

end


%%%%%%%%%%%%%%%%%%%%%%%%%% Corrector %%%%%%%%%%%%%%%%%%%%%%%%%%

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
% checkGradients_opts = optimoptions('fsolve',FiniteDifferenceType='forward'); checkGradients(Fcn,y,checkGradients_opts,Display='on',Tolerance=1e-6);   %FDM: Check the Jacobian matrix


% No initial solution found
if newton_flag~=1
    error_text = 'ERROR: No initial solution found!';                   % Set error text
    error_msg = 'CoSTAR stopped because corrector did not converge.';   % Set error message
    write_log(DYN,error_text)                                           % Write error text in log file
    write_log(DYN,'finalize_error',error_msg)                           % Finalize log file with error message
    if ~strcmpi(DYN.display,'off')
        disp(error_text); disp(' '); disp(error_msg); disp(' ')             % Display error text and message
        disp('-----------------------------------------------------')
        disp('--------------  Finished with error!  ---------------');
        disp('-----------------------------------------------------')
        disp(' ');
    end
    error(append(error_text(8:end),' ',error_msg))                      % Throw error
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Error Control Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
% This loop is slightly different from the loop in the Continuation class
% since it is does not update the dimension of the last curve point

added_output = cell(0,0);

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
            [y,~,newton_flag,~,J] = fsolve(Fcn,y0,newtonOpts);  % Solve corrector-function

        end

    end

    if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
        disp(' '); 
    end 
    write_log(DYN,'')
    added_output{1} = err;

end

info_text = 'Initial solution found!';


%%%%%%%%%%%%%%%%%%%%%%%%%% STABILITY %%%%%%%%%%%%%%%%%%%%%%%%%%
% Has to be called after the error_control step

if strcmpi(DYN.stability,'on')

    [multipliers,vectors,n_unstable,stability_flag] = ST.calc_stability(y,J);
    added_output{2} = multipliers;
    added_output{3} = n_unstable;
    added_output{4} = stability_flag;
    added_output{5} = vectors;

    if stability_flag == 0   
        warn_text = 'WARNING: Stability computation failed for the initial solution!';
        write_log(DYN,warn_text)
        S.warnings{end+1} = warn_text(10:end);                  % Save warning in Solution object
        warning(warn_text(10:end));
    elseif n_unstable == 0     
        info_text = 'Initial solution (stable) found!';
    elseif n_unstable > 0      
        info_text = 'Initial solution (unstable) found!';   
    end

end


%%%%%%%%%%%%%%%%%%%% Save initial solution %%%%%%%%%%%%%%%%%%%%

if ~strcmpi(DYN.display,'off')
    disp(info_text); disp(' '); 
end
write_log(DYN,info_text)

if strcmpi(DYN.sol_type,'quasiperiodic') && strcmpi(DYN.approx_method,'shooting')
    S.IF_arch_init_data(y,J,AM.Ik,newton_flag,AM.phi,DYN,added_output);
else
    S.IF_arch_init_data(y,J,newton_flag,DYN,AM,added_output);
end

% Check the frequencies
[warn_msg,~] = check_freq(DYN,y);
if ~isempty(warn_msg)                                                   % If frequency(s) are smaller than frequency limit
    DYN.cont = 'off';                                                   % Deactivate continuation (if it was set to 'on')
    warn_text = 'WARNING: Small or negative frequency(s) detected!';    % Set warning text
    S.stopping_flag = append(warn_msg);                                 % Save stopping message
    S.warnings{end+1} = warn_text(10:end);                              % Save warning in Solution object
    warning(warn_text(10:end));                                         % Display warning
    write_log(DYN,'finalize',append(warn_text,'\n\n',warn_msg))         % Finalize log file with warning text and message
    if ~strcmpi(DYN.display,'off')
        disp(' '); disp(warn_msg)                                       % Display warning text and message
        disp(' ')
        disp('-----------------------------------------------------')
        disp('------------  Finished with warning(s)!  ------------');
        disp('-----------------------------------------------------')
        disp(' '); disp(' ')
    end

% Set stopping message if no continuation is desired
elseif strcmpi(DYN.cont,'off')
    stopping_msg = 'CoSTAR stopped after initial solution because continuation was turned off.';    % Set stopping message
    write_log(DYN,'finalize',stopping_msg)                              % Finalize log file with stopping message
    S.stopping_flag = stopping_msg;                                     % Save stopping message
    if ~strcmpi(DYN.display,'off')
        disp(stopping_msg)                                                  % Display stopping message
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