% Function calculates initial curve point
%
% @DYN: DynamicalSystem object
% @S:   Solution object
% @AM:  ApproxMethod object
% @ST:  Stability object

function [S,AM,DYN] = initial_solution(DYN,S,AM,ST)

write_log(DYN,'initialize')                                             % Initialize log file

if ~strcmpi(DYN.display,'off')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%                                 CoSTAR                                  %')
    disp('%              Continuation of Solution Torus AppRoximations              %')
    disp(append('%                             Version ',DYN.costar_version,'                               %'))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' ')
end


%% Set function and options for fsolve

newtonOpts = optimoptions('fsolve','Display','iter-detailed','MaxFunEvals',1e5,'MaxIter',1e3,'FiniteDifferenceType','forward');

if strcmpi(DYN.approx_method,'shooting')
    newtonOpts.MaxIter = 50;
    newtonOpts.SpecifyObjectiveGradient = true;                         % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended
elseif strcmpi(DYN.approx_method,'finite-difference')                   % Special corrector function due to specification of Jacobian matrix
    newtonOpts.SpecifyObjectiveGradient = true;                         % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended    
elseif strcmpi(DYN.approx_method,'fourier-galerkin')
    % The following if-blocks must be executed here for the warning and info messages to be printed to the log (if they were executed in e.g. getIV, the log does not exist yet and the printing fails!)
    % Adapt n_FFT if not supplied, and if default value is not sufficient (if it is supplied by the user it was already checked in the gatekeeper)
    if ~isfield(DYN.opt_approx_method,'n_fft') && (max(AM.hmatrix,[],'all') >= AM.n_fft/2)
        n_fft_next_pow = nextpow2(2.*max(AM.hmatrix,[],'all')+1);  AM.n_fft = 2^n_fft_next_pow;
        info_text = sprintf('%s\n%s\n','At least one harmonic supplied in opt_approx_method.hmatrix violates the Nyquist-Shannon criterium in combination with the default value of n_fft = 2^6.', ...
                                        append('Therefore, the value of n_fft was adapted accordingly to 2^',num2str(n_fft_next_pow),' = ',num2str(AM.n_fft),'.'));
        write_log(DYN,info_text)
        if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'));  disp(info_text);  end
    end
    % If the error control is on and the maximum possible frequency (= n_fft/2-1, n_fft/2 is the Nyquist frequency) is already supplied, we need to increase AM.n_fft
    % Theoretically, we could continue without this adaption, but that would cause problems with the error control. See method "IF_increase_discretization" for details
    if strcmpi(AM.error_control,'on') && (max(abs(AM.hmatrix),[],'all') == (AM.n_fft/2-1))
        AM.n_fft = 2*AM.n_fft;
        info_text = sprintf('%s\n%s\n',append('The supplied hmatrix includes the highest possible harmonic (= n_fft/2-1 = ',num2str(0.5*AM.n_fft/2-1),') according to the Nyquist-Shannon criterium.'), ...
                            append('To avoid problems with the error control, the value of n_fft was increased from 2^',num2str(log2(0.5*AM.n_fft)),' to 2^',num2str(log2(AM.n_fft)),' = ',num2str(AM.n_fft),'.'));
        write_log(DYN,info_text)
        if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'));  disp(info_text);  end
    end
end

y0 = [AM.iv; DYN.param{DYN.act_param}];                                 % Initial value for the corrector
Fcn = @(y) AM.res_fun_init(y,y0);                                       % Function wrapper to set the complete residuum function for the initial solution


%% Corrector

warn_msg = cell(0,0);                                                   % Initialize

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


%% Check fsolve exit flags

error_flag = 1;                                                         % Initialise error_flag -> If set to 0: Critical error, CoSTAR is stopped
info_text = 'Initial solution found!';                                  % This message is updated right below or if stability is successfully computed

% No initial solution found or solution found but Jacobian can be undefined
if (newton_flag < 1) || (newton_flag == 2)
    error_flag = 0;                                                                                                         % Stop CoSTAR
    if newton_flag < 1
        error_text = append('ERROR: No initial solution found (fsolve exit_flag = ',num2str(newton_flag),')!');             % Set error text
        stopping_msg = 'CoSTAR stopped because corrector did not converge at initial solution.';                            % Set stopping message
    elseif newton_flag == 2
        error_text = 'ERROR: Equation solved, but change in y smaller is than the specified tolerance, or Jacobian at y is undefined (fsolve exit_flag = 2)!';     % Set error text
        stopping_msg = 'CoSTAR stopped because Jacobian can be undefined at initial solution.';                             % Set stopping message
    end

% Initial solution found, but fsolve exitflag is not default exitflag
elseif newton_flag == 3
    warn_msg{end+1} = 'WARNING: Equation solved for initial solution, but change in residual is smaller than specified tolerance (fsolve exit_flag = 3)!';
    S.warnings{end+1} = warn_msg{end}(10:end);                          % Save warning in Solution object
elseif newton_flag == 4
    warn_msg{end+1} = 'WARNING: Equation solved for initial solution, but magnitude of search direction is smaller than specified tolerance (fsolve exit_flag = 4)!';
    S.warnings{end+1} = warn_msg{end}(10:end);                          % Save warning in Solution object
end

% A user-defined Jacobian matrix can be checked here by executing the next line (since R2023b: function checkGradients is recommended instead of obj.fsolve_opts.CheckGradients = true)
% checkGradients_opts = optimoptions('fsolve',FiniteDifferenceType='forward'); checkGradients(Fcn,y,checkGradients_opts,Display='on',Tolerance=1e-5);   % Check the user-defined Jacobian matrix


%% Frequency Check

if error_flag == 1
    [stopping_msg,error_flag] = check_freq(DYN,y);
    if error_flag == 0                                                  % If frequency(s) are smaller than frequency limit
        error_text = 'ERROR: Small or negative frequency(s) detected for initial solution!';      % Set warning message
    end
end


    %% Error Control Loop
    % This loop is slightly different from the loop in the Continuation class
    % since it is does not update the dimension of the last curve point

    added_output = {NaN};                                               % Initialize (needed for saving)

    if strcmpi(AM.error_control,'on') && (error_flag == 1)

        iterate = 1;
        counter = 0;
        increase = 1;           % These values help to avoid a constant switching between de- and increasing, if error tolerance are set too close to each other
        decrease = 1;           % If we start increasing or decreasing once, the other one is not possible anymore

        while iterate

            counter = counter + 1;
            fieldnames_AM_struct = fieldnames(AM.ec_prop_save);                             % Dynamically save the properties of AM that are modified by the error control
            for i = 1:numel(fieldnames_AM_struct)                                           % In contrast to the error control in Continuation, we only have to save the AM properties here 
                AM.ec_prop_save.(fieldnames_AM_struct{i}) = AM.(fieldnames_AM_struct{i});   % We only need to reset some AM properties in the case that the discretization is decreased and the new error ...
            end                                                                             % is greater than the upper error limit (step width reduction is not possible here)
            err = AM.IF_estimate_error(y,DYN);
            err_control_text = append('Current Error: ',num2str(err),' -- Number of Harmonics: ',num2str(size(AM.hmatrix,2)-1));
            if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                disp(err_control_text);
            end
            write_log(DYN,err_control_text)

            %%%%%%%%%% Error is too high %%%%%%%%%%
            if (err > AM.error_limit(2)) && (counter < AM.ec_iter_max+1) && increase
                err_control_text = append('Increase Discretization (Iteration: ',num2str(counter),')');
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(err_control_text);
                end
                write_log(DYN,err_control_text)
                [y0,iterate] =  AM.IF_increase_discretization(y,DYN); %this might also update properties in AM
                decrease = 0;
            %%%%%%%%%% Error is too low %%%%%%%%%%
            elseif (err < AM.error_limit(1)) && (counter < AM.ec_iter_max+1) && decrease
                err_control_text = append('Decrease Discretization (Iteration: ',num2str(counter),')');
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(err_control_text);
                end
                write_log(DYN,err_control_text)
                [y0,iterate] = AM.IF_decrease_discretization(y,DYN); %this might also update properties in AM
                increase = 0;
            %%%%%%%%%% Error is fine %%%%%%%%%%
            else
                break           % Stop the loop since the error is fine
            end

            %%%%%%%%%% Recompute solution with in-/decreased discretization %%%%%%%%%%
            if iterate == 1
                % Corrector
                newtonOpts.Display = 'off';                                     % Deactivate fsolve output
                AM.IF_up_res_data(y0,DYN);                                      % Update AM properties and set y0 as initial value
                Fcn = @(y) AM.res_fun_init(y,y0);                               % Function wrapper to set the complete residuum function for the initial solution
                [y_ec,~,newton_flag_ec,~,J_ec] = fsolve(Fcn,y0,newtonOpts);     % Solve corrector-function

                % Exit flag handling Ia
                if any(newton_flag_ec == [1,3,4])       % Accept the result from fsolve if exitflag is 1, 3 or 4
                    % Check the new error: When decreasing the discretization, it can happen that the new error is above the upper error limit, which would require an increase again
                    % If this happens, the decreased solution y_ec is not accepted (since the error is too high). Instead: Reset the properties and exit the loop
                    if decrease
                        err_ec = AM.IF_estimate_error(y_ec,DYN);
                        if err_ec > AM.error_limit(2)                           % This nested if is done to avoid unnecessary calls of IF_estimate_error when increasing the discretization
                            for i = 1:numel(fieldnames_AM_struct)               % Reset the properties of AM that were modified by the error control -> needed to continue since new solution is not accepted
                                AM.(fieldnames_AM_struct{i}) = AM.ec_prop_save.(fieldnames_AM_struct{i});
                            end
                            ec_info_text = append('The current error ',num2str(err_ec),' is greater than the upper error limit of ',num2str(AM.error_limit(2)),['. ' ...
                                                  'Therefore, the new discretization (and the potential new n_fft value) is not accepted.']);
                            write_log(DYN,ec_info_text)
                            if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final')); disp(ec_info_text); end
                            decrease = 0;                                       % This line and the "continue" in the next line will a) print the old error and discretization again (for the sake of clarity) ...
                            continue                                            % and b) force the error control to go into the "Error is fine" block, which then terminates the error control loop
                        end
                    end
                    y = y_ec;
                    newton_flag = newton_flag_ec;
                    J = J_ec;
                    if newton_flag_ec == 3
                        warn_text = append('Equation solved for new discretization at initial solution, but change in residual is smaller than specified tolerance (fsolve exit_flag = 3)!');
                        write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
                        S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
                        warning(warn_text);                                                                         % Display warning
                    elseif newton_flag_ec == 4
                        warn_text = append('Equation solved for new discretization at initial solution, but magnitude of search direction is smaller than specified tolerance (fsolve exit_flag = 4)!');
                        write_log(DYN,append('WARNING: ',warn_text))                                                % Write warning in log file
                        S.warnings{end+1} = warn_text;                                                              % Save warning in Solution object
                        warning(warn_text);                                                                         % Display warning
                    end
                % Exit flag handling Ib
                else                                                % If exitflag from fsolve is < 0 or = 2
                    error_flag = 0;                                 % Stop CoSTAR
                    error_text = append('ERROR: Error control failed at initial solution for iteration ',num2str(counter),', because fsolve returned exit_flag = ',num2str(newton_flag_ec),'!');
                    stopping_msg = 'CoSTAR stopped because error control failed at initial solution.';  % Set stopping message
                    break                                           % Exit error control loop | Next: Display warning, write log file etc. -> done below
                end

            %%%%%%%%%% Increase of discretization failed %%%%%%%%%%
            elseif iterate == -1
                error_flag = 0;                                                    % Stop CoSTAR
                error_text = append('ERROR: Increase of discretization failed during error control at initial solution for iteration ',num2str(counter),'!');
                stopping_msg = 'CoSTAR stopped because error control failed at initial solution.';  % Set stopping message
                break                                                              % Exit error control loop | Next: Display warning, write log file etc. -> done below

            %%%%%%%%%% Increasing not possible: Maximum discretization reached %%%%%%%%%%
            elseif iterate == 0
                ec_info_text = 'Discretization cannot be increased any further since maximum discretization has been reached.';
                write_log(DYN,ec_info_text)                                        % Write info text in log file
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(ec_info_text);                                            % Display information
                end
                break                                                              % Exit error control loop (attention: error is ABOVE the upper limit)
            
            %%%%%%%%%% Decreasing not possible: Minimum discretization reached %%%%%%%%%%
            elseif iterate == 2
                ec_info_text = 'Discretization cannot be reduced any further since minimum discretization has been reached.';
                write_log(DYN,ec_info_text)                                        % Write info text in log file
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(ec_info_text);                                            % Display information
                end
                break                                                              % Exit the error control loop (error is below the bottom limit)

            %%%%%%%%%% Decrease of discretization failed %%%%%%%%%%
            elseif iterate == 3
                ec_info_text = append('Decrease of discretization failed. Continuing with previous discretization.');
                write_log(DYN,ec_info_text)                                        % Write info text in log file
                if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
                    disp(ec_info_text);                                            % Display information
                end
                break                                                              % Exit the error control loop (error is below the bottom limit)

            end

        end

        if ~(strcmpi(DYN.display,'off') || strcmpi(DYN.display,'final'))
            disp(' ');
        end
        write_log(DYN,'')
        added_output{1} = err;

    end


    %% Stability
    % Has to be called after the error_control step

    if strcmpi(DYN.stability,'on')

        if error_flag == 1
            [multipliers,vectors,n_unstable,stability_flag] = ST.calc_stability(y,J);
        else
            multipliers = NaN(DYN.dim,1);   vectors = NaN(DYN.dim);   n_unstable = NaN;   stability_flag = NaN;
        end
        added_output{2} = multipliers;
        added_output{3} = n_unstable;
        added_output{4} = stability_flag;
        added_output{5} = vectors;

        if stability_flag == 0                                                              % Stability computation failed
            warn_msg{end+1} = 'WARNING: Stability computation failed for initial solution!';% Set warning message
            S.warnings{end+1} = warn_msg{end}(10:end);                                      % Save warning in Solution object
        elseif n_unstable == 0                                                              % Stable solution found
            info_text = 'Initial solution (stable) found!';                                 % Update info text
        elseif n_unstable > 0                                                               % Unstable solution found
            info_text = 'Initial solution (unstable) found!';                               % Update info text
        end

    end


    %% Save initial solution 

    if strcmpi(DYN.sol_type,'quasiperiodic') && strcmpi(DYN.approx_method,'shooting')
        S.IF_arch_init_data(y,J,AM.Ik,newton_flag,AM.phi,DYN,added_output);
    else
        S.IF_arch_init_data(y,J,newton_flag,DYN,AM,added_output);
    end


%% Display info and warning and write log

if error_flag == 1                                  % Info (i.e. the message that a solution has been found): only if there is no critical error
    if ~strcmpi(DYN.display,'off')
        disp(info_text); disp(' ');                 % Display info text
    end
    write_log(DYN,info_text)                        % Write info text to log
end

if ~isempty(warn_msg)                               % Warning(s)
    for i = 1:numel(warn_msg)
        warning(warn_msg{i}(10:end)); disp(' ')     % Display warning
        if error_flag == 1
            write_log(DYN,append('\n',warn_msg{i})) % Write warning in log file
        else
            write_log(DYN,append(warn_msg{i},'\n')) % Write warning in log file
        end
    end
end


%% Stopping (no continuation)

if error_flag == 0                                  % Stopping in case of a critical error
    S.warnings{end+1} = error_text;                 % Save error text in Solution object
    write_log(DYN,error_text)                       % Write error text in log file
    write_log(DYN,'finalize_error',stopping_msg)    % Finalize log file with stopping message
    S.stopping_flag = stopping_msg;                 % Save stopping message
    warning(error_text); disp(' ');                 % Display error text
    if ~strcmpi(DYN.display,'off')
        disp(stopping_msg); disp(' ')               % Display stopping message
        disp('-----------------------------------------------------')
        disp('--------------  Finished with error!  ---------------')
        disp('-----------------------------------------------------')
        disp(' '); disp(' ')
    end
    DYN.cont = 'off';                               % Deactivate continuation (if it is enabled)

elseif strcmpi(DYN.cont,'off')                      % Stopping if there is no critical error
    stopping_msg = 'CoSTAR stopped after initial solution because continuation was turned off.';
    write_log(DYN,'finalize',stopping_msg)          % Finalize log file with stopping message
    S.stopping_flag = stopping_msg;                 % Save stopping message
    if ~strcmpi(DYN.display,'off')
        disp(stopping_msg); disp(' ')               % Display stopping message
        disp('-----------------------------------------------------')
        if ~isempty(lastwarn)
            disp('------------  Finished with warning(s)!  ------------')
        else
            disp('-------------- Successfully finished! ---------------')
        end
        disp('-----------------------------------------------------')
        disp(' '); disp(' ')
    end

end


%% Save DYN and S

if strcmpi(DYN.save,'on')
    save(append('CoSTAR_Save_',DYN.DYN_id,'.mat'),'DYN','S');        %save DYN and S in a mat-file
end


end