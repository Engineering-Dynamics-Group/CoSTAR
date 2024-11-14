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
        obj.p_error = AM.IF_estimate_error(obj.p_y1,DYN);
        err_control_text = append('Current Error: ',num2str(obj.p_error),' -- Number of Harmonics: ',num2str(size(AM.hmatrix,2)-1));
        write_log(DYN,err_control_text)
        if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
            disp(err_control_text); 
        end

        %Error is too high
        if (obj.p_error > AM.error_limit(2)) && (counter < AM.ec_iter_max+1) && increase
            [obj.yp,iterate] =  AM.IF_increase_discretization(obj.p_y1,DYN);                % This might also update properties in AM
            [obj.y0,obj.dy0] = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.y0,obj.dy0);      % Updates the older data points to match the dimension of the new solution vector (depends on solution type)
            obj.p_y0_old{2} = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.p_y0_old{2});
            obj.p_y0_old{3} = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.p_y0_old{3});
            decrease = 0;
            err_control_text = append('Increase Discretization (Iteration: ',num2str(counter),')');
            write_log(DYN,err_control_text)
            if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                disp(err_control_text); 
            end

        %Error is too low
        elseif (obj.p_error < AM.error_limit(1)) && (counter < AM.ec_iter_max+1) && decrease
            [obj.yp,iterate]  = AM.IF_decrease_discretization(obj.p_y1,DYN);                % This might also update properties in AM 
            [obj.y0,obj.dy0]  = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.y0,obj.dy0);     % Updates the older data points to match the dimension of the new solution vector (depends on solution type)
            obj.p_y0_old{2} = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.p_y0_old{2});
            obj.p_y0_old{3} = AM.IF_update_sol_dim(DYN,numel(obj.yp),obj.p_y0_old{3});
            increase = 0;
            err_control_text = append('Decrease Discretization (Iteration: ',num2str(counter),')');
            write_log(DYN,err_control_text)
            if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                disp(err_control_text); 
            end 

        else
            iterate = 0;
        end

        % Recalculate a new curve point with de-/increased discretisation
        if iterate  
            if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
                obj.fsolve_opts.MaxIter = 30;
                AM.IF_up_res_data(obj,DYN);
                AM.y_old = obj.yp;                                                          % If autonomous y_old contains predicted solution
                Fcn = @(y)AM.fun_Jac_wrapper(y,DYN,obj);                                    % Set functionwrapper to provide Jacobian
                obj.fsolve_opts.SpecifyObjectiveGradient = true;
            else
                AM.IF_up_res_data(obj);                                                     % Archive initial solution
                Fcn = @(y)[AM.res(y);obj.sub_con(y,obj)];                                   % Define corrector-function containing the residual function and the subspace-constraint
            end

            [obj.p_y1,~,obj.p_newton_flag,obj.p_output,obj.p_J1] = fsolve(Fcn,obj.yp,obj.fsolve_opts);  % Solve corrector-function

        end

    end
    
  
end