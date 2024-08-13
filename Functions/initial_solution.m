% Function calculates initial curve point to start the continuation toST
% given residual-equation
%
% @DYN: DynamicalSystem object
% @S:   Solution object
% @AM:  ApproxMethod object
% @ST:  Stability object

function [S,AM,DYN] = initial_solution(DYN,S,AM,ST)

newtonOpts = optimoptions('fsolve','Display','iter-detailed','MaxFunEvals',1e5,'FiniteDifferenceType','forward');
y0 = [AM.iv;DYN.param{DYN.act_param}];
if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
    Fcn = @(y)AM.fun_Jac_wrapper_init(y,y0,DYN);                    % Function wrapper for initial solution, if Jacobian is supplied
    newtonOpts.SpecifyObjectiveGradient = true;
elseif strcmpi(DYN.approx_method,'finite-difference')               % Special corrector function due to specification of Jacobian matrix
    newtonOpts.SpecifyObjectiveGradient = true;                     % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended
    Fcn = @(y) AM.corr_fun_init_FDM(y,y0);
else
    Fcn = @(y)[AM.res(y);y(end)-y0(end)];                           % Last entry is not necessary, but this way the Jacobian has the correct dimension
end


%%%%%%%%%%%%%%%%%%%%%%%%%% Corrector %%%%%%%%%%%%%%%%%%%%%%%%%%

[y,~,newton_flag,~,J] = fsolve(Fcn,y0,newtonOpts);
% checkGradients_opts = optimoptions('fsolve',FiniteDifferenceType='forward'); checkGradients(Fcn,y,checkGradients_opts,Display='on',Tolerance=1e-6);   %FDM: Check the Jacobian matrix

if newton_flag~=1                           % If no initial solution is beeing found, throw error message
    S.stopping_flag = 'CoSTAR stopped because corrector did not converge.';
    error('No initial solution found!')
else
    S.stopping_flag = 'CoSTAR stopped after initial solution because continuation was turned off.';
    [warn_msg,~] = check_freq(DYN,y);       % Check the frequencies
    if ~isempty(warn_msg)                   % If frequency(s) are smaller than frequency limit
        S.stopping_flag = warn_msg;
        DYN.cont = 'off';                   % Deactivate continuation (if it was set to 'on')
        warning(warn_msg)                   % Throw warning message
    else


        %%%%%%%%%%%%%%%%%%%%%%%%%% Error Control Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
        %This loop is slightly different from the loop in the Continuation class
        %since it is does not update the dimension of the last curve point

        added_output = cell(0,0);

        if strcmpi(AM.error_control,'on')

            iterate = 1;
            counter = 0;
            increase = 1;           % These values help to avoid a constant switching between de- and increasing, if error tolerance are set too close to each other
            decrease = 1;           % If we start increasing or decreasing once, the other one is not possible anymore

            while iterate

                counter = counter +1;
                err = AM.IF_estimate_error(y,DYN);
                disp(append('Current error: ',num2str(err)));

                if (err > AM.error_limit(2)) && (counter < AM.ec_iter_max+1) && increase
                    disp(append('Increase Discretization. Iteration ',num2str(counter)));
                    [y0,iterate] =  AM.IF_increase_discretization(y,DYN); %this might also update properties in AM
                    decrease = 0;

                    %             elseif (err < AM.error_limit(1)) && (counter < AM.ec_iter_max+1) && decrease
                    %
                    %                 [y0,iterate] = AM.IF_decrease_discretization(y,DYN); %this might also update properties in AM
                    %                 increase = 0;
                else
                    iterate = 0;
                end

                if iterate  %Recalculate a new curve point with de-/increased discretization
                    %If autonomous y_old contains predicted solution
                    if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
                        AM.IF_up_res_data(obj,DYN); %!!! Check this line if error_control becomes available for quasiperiodic shooting !!!
                        AM.y_old = y0;
                        Fcn = @(y)AM.fun_Jac_wrapper_init(y,y0,DYN);                                            %Function wrapper for initial solution, if Jacobian is supplied
                        newtonOpts.SpecifyObjectiveGradient = true;
                    else

                        AM.iv = y0(1:(end-1));                                                         %Archive initial solution
                        Fcn = @(y)[AM.res(y);y(end)-y0(end)];                                           %define corrector-function containing the residual function and the subspace-constraint
                    end

                    [y,~,newton_flag,~,J] = fsolve(Fcn,y0,newtonOpts);  %solve corrector-function

                end
            end

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
        end


        %%%%%%%%%%%%%%%%%%%% Save initial solution %%%%%%%%%%%%%%%%%%%%

        disp('-------------------------------------------------------')
        disp('--------------- Initial solution found! ---------------')
        disp('-------------------------------------------------------')

        if strcmpi(DYN.sol_type,'quasiperiodic') && (strcmpi(DYN.approx_method,'shooting') || strcmpi(DYN.approx_method,'mshm'))
            S.IF_arch_init_data(y,J,AM.Ik,newton_flag,AM.phi,DYN,added_output);
        else
            S.IF_arch_init_data(y,J,newton_flag,DYN,AM,added_output);
        end

    end

end

end