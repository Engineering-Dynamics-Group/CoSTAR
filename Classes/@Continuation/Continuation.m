% Class continuation provides a path continuation algorithm to continue
% curves defined by the roots of an algebraic system of equations
classdef Continuation < handle

    properties
        cont = 'on';                                                        %Continue if cont='on', if cont='off' stop after calculation of initial solution (actually not used here)
        pred string = 'tangent';                                            %Predictor to be used for the prediction of new curve point
        subspace string = 'pseudo-arc';                                     %Defines the subspace-constraint to close the corrector-equation
        direction = 1;                                                      %Direction in which the curve shall be continued. 1 = positive mu-direction, -1 = negative mu-direction
        plot = 'on';                                                        %If = 'on', a continuation plot is displayed during continuation
        step_control = 'angle';                                             %Default method used by step control
        step_control_param;                                                 %Array storing parameters used by step control
        step_width = 0.1;                                                   %Step width for continuation, further the initial step width if step control is implemented
        step_width_limit = [];                                              %Limits for adaption of step width by step control
        max_cont_step = 1500;                                               %Maximum number of continuation steps
        mu_limit                                                            %Limit of bifurcation parameter, when reached the continuation stops, e.g. mu_limit = [-1,2]
        fsolve_opts = optimoptions('fsolve','Display','none','UseParallel',false,'MaxFunctionEvaluations',20000,'MaxIter',1e3); %Options for fsolve to solve corrector-equation
        yp                                                                  %YP is the predicted curve point. It is public since it is exchanged with the ApproxMethod object.
        sub_con function_handle                                             %Function_handle for subspace constraint 
        d_sub_con function_handle                                           %Derivative of subspace constraint   

        %Parameter of active curve point
        dy0                                                                 %Direction vector of the predictor
        y0                                                                  %Curve point

    end
    
    properties(SetAccess = private, GetAccess = public)                     %These properties are stored into the solution objects
        

        %Parameters of predicted curve point
        p_y1 = 0;                                                             %Initialise for get method of p_arcl_1 to be well defined
        p_J1
        p_newton_flag = 0;                                                    %Exit flag of Newton solver
        p_stopping_flag                                                       %Exit flag of continuation
        p_arcl_0  = 0;                                                        %arc-length of current point
        p_arcl_1  = 0;                                                        %arc-length of new point: Has its own get method

        p_y0_old                                                              %Stores last 3 "old" curve point for predictors

        p_error = NaN;                                                        %Current error of the solution point
        p_n_unstable_1 = NaN;                                                 %Indicates the number of unstable multipliers
        p_n_unstable_0 = NaN;                                                 %Indicates the number of unstable multipliers curve point
        p_multipliers = NaN;                                                  %Multipliers to determine stability (e.g. Eigenvalues, Floquet-Multipliers or Lyapunov-Exponents)
        p_max_multiplier_1 = NaN;                                             %Maximal multiplier responsible for stability loss or gain of current curve point
        p_max_multiplier_0 = NaN;                                             %Maximal multiplier responsible for stability loss or gain of last curve point
        p_vectors                                                             %(Eigen-)vectors corresponding to multipliers

        %Parameters of bifurcation points
        p_y_bfp                                                               %curve point at bifurcation point
        p_J_bfp                                                               %Jacobian matrix at bifurcation point
        p_multipliers_bfp                                                     %multipliers at bifurcation point
        p_vectors_bfp                                                         %(Eigen-)vectors corresponding to multipliers at bifurcation point
        p_arclength_bfp                                                       %arc-length of current bifurcation point 
        p_error_bfp  = NaN;                                                   %error at the bifurcation point 
        p_stability_flag                                                      %exitflag of stability computation
        p_newton_flag_bfp                                                     %fsolve exitflag of iterating bifurcation point
    end
    
    properties(Access=private)
        p_contDo uint32 = 1;                                                  %While contDo=1, do continuation, if conDo=0 stop continuation
        p_local_cont_counter uint32 = 1;                                      %Counts number of computed solutions
        p_last_msg;                                                           %Saves the latest "Iter" messages printed to the command window

        %Parameters of the last point
        p_dy_old                                                              %Direction vector of the predictor
        p_r_old = 1;                                                          %Factor which adapts step width
        p_dx_dmu_old                                                          %Used for PID step control
        p_e_old = 1;                                                          %Used for PID step control
        p_axes_values_old                                                     %Axes values of continuation plot of old solution
        p_stability_flag_old                                                  %Old exitflag of stability computation
                
        %Parameters of penultimate point
        p_e_old_old = 1;                                                      %Used for PID step control

        %Parameter of active curve point
        p_mu0                                                                 %Continuation parameter
        p_J0                                                                  %Jacobian
        p_period0%%%                                                          %Period of Solution
        p_it                                                                  %Number of iterations
        p_output                                                              %Output information of Newton-solver (fsolve)
        p_convergence = 1;                                                    %Set to zero if fsolve did not converge and step_width is reduced
        p_step_width_init                                                     %Storage for initial step_width
        p_limit                                                               %needed for plotting
        p_r = 1;                                                              %Factor which adapts step width
        p_e = 1;                                                              %Used for PID step control
        p_dx_dmu                                                              %Used for PID step control

        p_initial_slope                                                       %Initial slope for secant predictor (only used in initial continuation step)
        p_use_qr logical = false;                                             %Boolean to determine whether a qr decomposition has to be used to determine inital slope (only if secant method failed)

    end
    %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%
    
    methods
        %% Constructor
        function obj = Continuation(options)    

            obj = updateoptions(obj,options);
            if strcmpi(obj.step_control,'on')                               %If user sets options.opt_cont.step_control = 'on', ...
                obj.step_control = 'angle';                                 %the default step control method 'angle' is chosen
            end 
            obj = choose_stepcontrol_param(obj);                            %Set default values for step_control_param
            obj = choose_subspace(obj);                                     %Define subspace_constraint and derivative

            obj.p_step_width_init = obj.step_width;                         %Set initial step_width
            %Set the step width limits if they were not given by user
            if isempty(obj.step_width_limit); obj.step_width_limit = [0.2.*obj.step_width,5.*obj.step_width]; end
        
        end

        %% Main method
        S = m_continuation(obj,DYN,S,AM,ST);                                %Main method of class Continuation
        
        %% Methods
        obj = initial_slope(obj,DYN,AM);                                    %Calculates second curve point to determine initial slope for direction vector (secant, parable, cubic)
        obj = direction_vector(obj);                                        %Calculates the vector of the direction of the predictor (tangent, secant, ...)
        obj = predictor(obj);                                               %Generates predictor to given point on diagram
        
        obj = iterate_data(obj);                                            %Iterates the relevant data: new data point 1 is now current data point 0 for next iteration
        obj = check_limits(obj,DYN);                                        %Checks, if the current continuation parameter is within the prescribed mu limit

        obj = stepcontrol(obj,DYN);                                         %Method adapts step width according to achieve a determined number of Newton iterations
        obj = choose_stepcontrol_param(obj);                                %Method defining the default values of step_control_param
        obj = choose_subspace(obj);                                         %Method defines function for subspace constraint
        
        obj = plot_contplot(obj,S,DYN);                                     %Method for plotting a live continuation plot
        obj = error_control(obj,S,AM,DYN);                                  %Method for controlling the error by adapting the discretisation scheme
        obj = bifurcation_stability(obj,DYN,AM,S,ST);                       %Method for calculating the stability and determining bifurcation points

        %% Get Method
        %By altering the get-method. p_arcl_1 always computes the current arc-length
        function p_arcl_1 = get.p_arcl_1(obj)
            p_arcl_1 = obj.p_arcl_0 + obj.direction.*norm(obj.p_y1-obj.y0); %This is the current arc length
        end


    end
    %%%%%%%%%%%%%%%%
    methods(Static)
        s_CON_gatekeeper(GC,opt_con);                                       %Gatekeeper for the opt_con option structure
        help_struct = s_help_opt_cont();                                    %Help function supplies the help struct for the opt cont option structure. 
    end
    
end