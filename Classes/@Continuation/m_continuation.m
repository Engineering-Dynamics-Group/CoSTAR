% Main method which continues a curve defined by the roots of an algebraic
% system of equations. This function needs an object of Solution (S) and an
% object of ApproxMethod (AM). m_continuation returns an object of Solution
% which contains the continued curve
%
%@obj:      Continuation class object
%@DYN:      DynamicalSystem object
%@S:        Solution subclass object
%@AM:       ApproxMethod sublass object
%@ST:       Stability subclass object
%
%@S:        Solution subclass object

function  S = m_continuation(obj,DYN,S,AM,ST)


%%%%%%%%%%%%%%%  INITIALISATION  %%%%%%%%%%%%%
obj.y0  = S.y0;                                                         %Get active curve point from object of Solution class
obj.p_y0_old = num2cell([zeros(numel(obj.y0),2),obj.y0],1);             %Must be initialized as 1x3 cell array for the error control to work
obj.p_mu0 = S.y0(end,1);
if isa(S.J,'cell'); obj.p_J0  = S.J{1,1}; else; obj.p_J0  = S.J; end
if strcmpi(DYN.stability,'on') && strcmpi(ST.iterate_bfp,'on')
    obj.p_n_unstable_0 = S.n_unstable; 
    ST.update_curve_container(DYN,AM,S.arclength,obj.y0,S.multipliers,obj.p_n_unstable_0); %Update the curve container with the first point
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  CONTINUATION  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while  obj.p_contDo
        
    %%%%%%%%%%%%%%%%%%  PREDICTOR AND STEP CONTROL  %%%%%%%%%%%%%%%%%
    obj = obj.direction_vector();                   %calculate direction vector
   
    if obj.p_newton_flag >= 1                       %stepcontrol may be called if corrector converged (first loop: obj.p_newton_flag=0 as stepcontrol does not make sense)
        obj = obj.stepcontrol();                    %adapt step width
        obj.p_convergence = 1;                      %reset convergence property if corrector did not converge previously
    end     
    
    obj = obj.predictor();                          %calcuate predicted point

    [warn_msg,~] = check_freq(DYN,obj.yp);          %check the frequencies at the predictor point (this is not done in method predictor in order to be able to break the while loop)
    if ~isempty(warn_msg)                           %if frequency(s) are smaller than frequency limit
        obj.p_stopping_flag = 3;                    %set termination flag
        warning(warn_msg)                           %throw warning message
        break                                       %break the continuation loop (any further computation below might fail)
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%  CORRECTOR  %%%%%%%%%%%%%%%%%%%%%%%%%
    if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))                %special corrector function for quasi-periodic shooting
        obj.fsolve_opts.MaxIter = 50;
        AM.IF_up_res_data(obj,DYN);                                                                 %pass information to the ApproxMethod object
        AM.y_old = obj.yp;                                                                          %if autonomous y_old contains predicted solution
        Fcn = @(y)AM.fun_Jac_wrapper(y,DYN,obj);                                                    %set functionwrapper to provide jacobian
        obj.fsolve_opts.SpecifyObjectiveGradient = true;                                            %Jacobian matrix is passed by the user
    elseif strcmpi(DYN.approx_method,'finite-difference')                                           %special corrector function for FDM due to specification of Jacobian matrix
        AM.IF_up_res_data(obj);                                                                     %pass information to the ApproxMethod object
        obj.fsolve_opts.SpecifyObjectiveGradient = true;                                            %Jacobian matrix is passed by the user | obj.fsolve_opts.CheckGradients = true; can be used to automatically check the Jacobian matrix
        Fcn = @(y) AM.corr_fun_FDM(y,obj);                                                          %set corrector-function
    else
        AM.IF_up_res_data(obj);                                                                     %pass information to the ApproxMethod object
        Fcn = @(y)[AM.res(y);obj.sub_con(y,obj)];                                                   %define corrector-function containing the residual function and the subspace-constraint
    end
    
    [obj.p_y1,~,obj.p_newton_flag,obj.p_output,obj.p_J1] = fsolve(Fcn,obj.yp,obj.fsolve_opts);      %solve corrector function


    %%%%%%%%%%%%%%%%%%  IF NO CONVERGENCE OF FSOLVE  %%%%%%%%%%%%%%%%
    if(obj.p_newton_flag<1 && obj.step_width<=obj.step_width_limit(1,1))                            %if fsolve did not converge and step width is already <= minimal step width 
        obj.p_stopping_flag = 0;                                                                    %set exitflag
        obj.p_contDo = 0;                                                                           %stop continuation
        warning('Continuation stopped because corrector did not converge and step width has reached minimal value!');
    
    elseif(obj.p_newton_flag<1 && obj.step_width>obj.step_width_limit(1,1))                         %if fsolve did not converge and step width is above minimal step width 
        step_width_pre = 0.5.*obj.step_width;                                                       %new preliminary step width
        obj.step_width = max([step_width_pre,obj.step_width_limit(1)]);                             %set step_width. If new preliminary step width falls below minimal step width, take minimal step width
        obj.p_convergence = 0;                                                                      %set property p_convergence to zero (for resetting the step_width after convergence)
        disp(append('Stepwidth adapted to stepwidth = ',num2str(obj.step_width),', because corrector did not converge!'));      
    

    else
        %FDM: Check the Jacobian matrix -> Since R2023b: checkGradients is recommended instead of obj.fsolve_opts.CheckGradients = true
        % checkGradients_opts = optimoptions('fsolve',FiniteDifferenceType='forward'); checkGradients(Fcn,obj.p_y1,checkGradients_opts,Display='on',Tolerance=1e-6);

        %%%%%%%%%%%%%%%%%%%%%%%%%  FREQUENCY CHECK  %%%%%%%%%%%%%%%%%%%%%%%
        [warn_msg,~] = check_freq(DYN,obj.p_y1);    %check the frequencies of the solution
        if ~isempty(warn_msg)                       %if frequency(s) are smaller than frequency limit
            obj.p_stopping_flag = 3;                %set termination flag
            obj.p_contDo = 0;                       %stop continuation
            warning(warn_msg);                      %throw warning message
            
        else

            %%%%%%%%%%%%%%%%%%%%%%%%%  ERROR CONTROL  %%%%%%%%%%%%%%%%%%%%%%%
            %Control the error by adapting the discretization. If the ansatz function is adapted, error_control solves the equation system again
            %and computes new values, which get saved and iterated, which is why the error_control must be called before IF_arch_data and iterate_data.
            if strcmpi(AM.error_control,'on'); obj.error_control(S,AM,DYN); end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%  STABILITY  %%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate the Lyapunov stability. If a bifurcation point is found, the location of the point can be iterated
            %IMPORTANT: Has to be called after the error_control step
            if strcmpi(DYN.stability,'on'); obj.bifurcation_stability(DYN,AM,S,ST); end


            %%%%%%  SAVING, PLOTTING, CHECKING LIMITS AND DISPLAY INFO  %%%%%
            %IMPORTANT: order of calling the methods must not be changed
            S.IF_arch_data(obj,DYN,AM);                                             %store calculated data point in Solution object

            if strcmpi(obj.plot,'on'); obj.plot_contplot(S,DYN); end                %display a continuation plot if desired

            obj = obj.iterate_data();                                               %store calculated data point as current data point for next continuation step

            obj = obj.check_limits();                                               %check limits for bifurcation parameter and number of steps and display information

        end

    end 

end

switch obj.p_stopping_flag                      %Set stopping flag of continuation
    case 0
        S.stopping_flag = 'CoSTAR stopped because corrector did not converge and step width has reached minimal value.';
    case 1
        S.stopping_flag = 'CoSTAR stopped because limit of continuation parameter was reached.';
    case 2
        S.stopping_flag = append('CoSTAR stopped because maximal number of continuation steps (', num2str(obj.max_cont_step), ') was reached.');
    case 3
        S.stopping_flag = warn_msg;
end

end