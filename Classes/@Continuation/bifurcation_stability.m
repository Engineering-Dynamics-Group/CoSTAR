%This function determines the stability by calling the appropriate function from the Stability subclass
%and determines bifurcations points. If bifurcation points are present, they are iterated 
%
%@obj:  Continuation class object
%@DYN:  Dynamical system class object
%@AM:   Approximation Method subclass object
%@S:    Solution subclass object
%@ST:   Stability subclass object

function bifurcation_stability(obj,DYN,AM,S,ST)


        %% Calculate Stability
        [obj.p_multipliers,obj.p_vectors,obj.p_n_unstable_1,stability_flag] = ST.calc_stability(obj.p_y1,obj.p_J1);               % Compute the respective multiplier (eigenvalue, Floquet, Lyapunov Exponent)
        

        if stability_flag > 0                                           % Only move forward if the stability computation was regular
        

            %% Iteration of bifurcation point
            if strcmpi(ST.iterate_bfp,'on')
        
                ST.update_curve_container(DYN,AM,obj.p_arcl_1,obj.p_y1,obj.p_multipliers,obj.p_n_unstable_1);       % Save the current solution curve point into the container (FIFO storage of the last n points)
        
                if ~(obj.p_n_unstable_0==obj.p_n_unstable_1)            % Change in number of unstable multipliers detected
        
                    iterate = 1;
                    counter = 0;
                    multipliers_old = obj.p_multipliers;                % For Evaluating stopping criteria
                    

                    while iterate

                        counter = counter + 1;
        
                        % Initial point
                        obj.yp = ST.approx_posc(DYN);                   % Function approximates the point of stability change based on the curve_container; output needs to be assigned to yp for the up_res_data method to work
        
                        % Correct
                        if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
                            AM.IF_up_res_data(obj,DYN);
                            Fcn = @(y)AM.fun_Jac_wrapper(y,DYN,obj);                % Set functionwrapper to provide Jacobian
                            obj.fsolve_opts.SpecifyObjectiveGradient = true;
                        elseif strcmpi(DYN.approx_method,'finite-difference')
                            AM.IF_up_res_data(obj);
                            obj.fsolve_opts.SpecifyObjectiveGradient = true;
                            Fcn = @(y) AM.corr_fun_FDM(y,obj);
                        else
                            AM.IF_up_res_data(obj);                                 % Archive initial solution
                            Fcn = @(y)[AM.res(y);obj.sub_con(y,obj)];               % Define corrector-function containing the residual function and the subspace-constraint
                        end
                        [obj.p_y_bfp,~,p_newton_flag_bfp,~,obj.p_J_bfp] = fsolve(Fcn,obj.yp,obj.fsolve_opts);               % Solve corrector-function
        
                        % Get Stability
                        [obj.p_multipliers_bfp,obj.p_vectors_bfp,n_unstable,stability_flag] = ST.calc_stability(obj.p_y_bfp,obj.p_J_bfp);     % Compute the stability of the newly found point

                        % Update curve_container
                        ST.update_curve_container_bfp(obj.p_y_bfp,obj.p_multipliers_bfp,n_unstable);                        % Update the curve container data with the newly found point at the appropriate position
        
                        % Stopping criteria
                        if (max(abs(obj.p_multipliers_bfp-multipliers_old))<ST.abstol_multiplier) || (counter>(ST.max_iter-1)) || p_newton_flag_bfp<1 || stability_flag<1
                            iterate = 0;
                        end
                        
                        multipliers_old = obj.p_multipliers_bfp;                    % For Evaluating stopping criteria
        
                    end

                    
                    % Save the bifurcation point
                    if (p_newton_flag_bfp>0) && (stability_flag>0)
                        
                        % Get the arc-length of the bifurcation point (for solution object)
                        [~,idx] = min(abs(cell2mat(ST.curve_container(3,:))));      % This should identify the (approximated) bifurcation point
                        obj.p_arclength_bfp = ST.curve_container{2,idx};            % This is the arclength of the (approximated) bifurcation point
                        
                        if strcmpi(AM.error_control,'on')
                            obj.p_error_bfp = AM.IF_estimate_error(obj.p_y_bfp,DYN); 
                        end

                        % IF_arch_bfp_data MUST be called before clean_curve_container for the identify_bifurcation method to work properly
                        S.IF_arch_bfp_data(obj,DYN,AM,ST);                          % Stores the data of the iterated bifurcation point
                    
                    else

                        warning('Iteration of the bifurcation point failed on the current step!');

                    end

                    ST.clean_curve_container;                                       % Deletes all elements before the stability change
        
                end
        
            end


        %% Warning if stability computation has failed
        else

            warning('Stability computation failed on the current step!');

        end


end