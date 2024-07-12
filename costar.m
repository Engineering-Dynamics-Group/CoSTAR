%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% This function continues stationary solutions of dynamical systems. These
% can be either stable-state solutions or periodic solutions.  
%
% Input Arguments:
%
% @options:     Struct containing all options for calculating solution. Datatype of all fields of @options have to be structs itself
%
%               Mandatory fields (datatype: struct):
%
%                       'system':           Mandatory fields of 'system':
%                                               --> 'order':                Order of ODE (0 if algebraic equation) 
%                                               --> 'rhs':                  Right hand side of ode or ae
%                                               --> 'dim':                  Dimension of the state space                                                                          
%                                           Optional fields of 'system':
%                                               --> 'param':                Vector of parameters of the system 
%                                               --> 'info':                 Information about the system (string)
%
%                       'opt_sol':          Mandatory fields of 'opt_sol':
%                                               --> 'sol_type':             Type of solution to be continued ('equilibrium','periodic')  
%                                               --> 'cont':                 Option whether solution shall be continued (1: continuation, 0: single solution)
%                                           Optional fields of 'opt_sol':
%                                               --> 'approx_method':        ApproxMethod specific method to calculate solution (e.g. for sol_type = 'periodic' => approx_method = 'shoot')
%                                               --> 'act_param':            Parameter which shall be continued (e.g. param = [lambda,mu,eta], mu is bifurcation parameter => act_param = 2)
%                                               --> 'non_auto_freq':        Value or function_handle as function of bifurcation parameter for that returns the non-autonomous frequency (e.g. @(mu) mu)
%                                               --> 'auto_freq':            If system is autonomous set auto_freq to 1
%
%                       'opt_init':         Mandatory fields of 'opt_init':
%                                               --> 'ic':                   Initial condition for continuation
%                                           Optional fields of 'opt_init':
%                                               -- No optional fields existing -- 
%
%               Optional fields (datatype: struct):
%
%                       'opt_approx_method':                                Depends on chosen solution method
%
%                       'opt_cont':         Mandatory fields of 'opt_cont':
%                                               -- No mandatory fields existing -- 
%                                           Optional fields of 'opt_cont':
%                                               --> 'pred':                 Method to calculate predictor point (default: 'tangent'; allowed: 'tangent', 'secant') [secant might be faster for larger systems] 
%                                               --> 'subspace':             Solution subspace constraint (default: 'pseudo-arc'; allowed: 'tangent','arclength','natural')
%                                               --> 'direction':            Direction of continuation (value: 1 or -1) [to be replaced in further version]
%                                               --> 'stability':            No stability method implemented yet
%                                               --> 'maxcontstep':          Maximum number of continuation steps (default: 1500)
%                                               --> 'step_width':           Step width for continuation in case of no step control. If SC is enabled, it sets the initial step width (default: 0.1)
%                                               --> 'stepcontrol':          Method used by step control (default: 'angle'; allowed: 'off', 0, 'on', 1, 'corrector_iterations', 'norm_corrector_predictor', 'angle', 'combination', 'pid')
%                                               --> 'step_width_limit'      Boundaries for step with when step control is used (default: [0.2, 5].*initial_step_width)
%                                               --> 'it_nominal'            Number of corrector iterations which shall be achieved using step control (default: 3)
%                                               --> 'plot':                 Activates the live plot of the continued curve (value 1 or 0; default: 1)
%                                               --> 'mu_limit':             Boundaries for continuation (e.g. [1,5])
%                                               --> 'display':              Whether a console output during the continuation is to be made (values: 0 or 1)
%
%
% Output Arguments: 
%
% @S:           Contains all simulated data
%                                           fields:
%                                               --> 'freq':                 List of frequency of solution at curve point
%                                               --> 'y0':                   Initial curve point
%                                               --> 's':                    Solution vector to curve point
%                                               --> 'J':                    Jacobian matrix to curve point
%                                               --> 'mu':                   List of continuation parameters to curve points
%                                               --> 'dy':                   Tangent vector of curve point
%                                               --> 'multipliers':
%                                               --> 'newton_flag':
%                                               --> 'flag':                 Flag for continutation algorithm (0 = continuation not successful, 1 = continuation successful, 2 = continuation terminated by number of curve points)
%                                               --> 'step_width':           Step width from curve point to next curve point
%                                               --> 'arclength':            Arclength of curve upt to curve point
%                                               --> 'stability':            Logical values for stability of solution at continued curve point
%                                               --> 'S_id':                 Identifier
%
% @DYN:         Contains all information about the system
%                                           fields:
%                                               --> 'info':                 Information about the system
%                                               --> 'order':                Order of ODE (0 if algebraic equation)
%                                               --> 'rhs':                  Right-hand-side of ODE or AE
%                                               --> 'param':                Vector of parameters
%                                               --> 'sol_type':             SolutionType (e.g. periodic)
%                                               --> 'approx_method':        Solution method (e.g. Shoot)
%                                               --> 'cont':                 Use path continuation (=1)
%                                               --> 'act_param':            Active parameter of parameter vector
%                                               --> 'non_auto_freq':        Non-autonomous frequencies
%                                               --> 'auto_freq':            Autonomous frequencies (initial values)
%                                               --> 'ic':                   Initial condition
%                                               --> 'iv':                   Initial value
%                                               --> 'system':               systems struct
%                                               --> 'opt_sol':              Solution struct
%                                               --> 'opt_approx_method':
%                                               --> 'opt_init':
%                                               --> 'opt_cont':
%                                               --> 'n':
%                                               --> 'n_auto':
%                                               --> 'n_freq':
%                                               --> 'DYN_id':
%
%
% Example:   
%               param = [0.1,0.3,0.5];
%               Fcn = @(t,z,param)duffing_ap(t,z,param);
%
%               options.system   = costaropts('order',1,'rhs',Fcn,'param',param);    
%               options.opt_sol  = costaropts('cont',1,'non_auto_freq',@(mu) mu,'sol_type','periodic','approx_method','shooting','act_param',2); 
%               options.opt_init = costaropts('ic',IC);
%               options.opt_cont = costaropts('pred','tangent','subspace','pseudo-arc','mu_limit',[1,3.5],'stepcontrol',1);             
%               options.opt_approx_method = costaropts('solver','ode45'); 
%           
%               [S,DYN] = costar(options);
%
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function [S,DYN] = costar(options)

    %% Gatekeeper
    GC = Gatekeeper();                                                         %Create the gatekeeper.
    options = GC.m_gatekeeper(options);                                        %This main function initially controls all input data by calling various static methods of the classes DynamicalSystem, SolutionType and Cotinuation    

    clear GC;                                                                  %Gatekeeper has done it's job.
    %% Dynamical System
    DYN = DynamicalSystem(options);

    %% Approximation Methods
    AM = ApproxMethod.s_method_selection(DYN);                                  %Calls the static method_selection methods of the SuperClass SolutionType, which creates the appropriate object of a subclass (Shoot, etc.)
    
    %% Create object of Solution
    S = Solution.s_solution_selection(DYN,AM);                                  %Calls the static solution_selection method of the SuperClass Solution, which creates the appropriate object of a subclass (EQ_Sol, PS_Shoot_Sol etc.)
    
    %% Stability class                                                          
    ST = Stability.s_stability_selection(DYN,AM);                           %Calls the static method_selection method of the SuperClass Stability,   which creates the approriate object of a subclass
    %% Calculate initial solution
    [S,AM,DYN] = initial_solution(DYN,S,AM,ST);

    %% Continuation
    if strcmpi(DYN.cont,'on')
        CON = Continuation(options.opt_cont);
        S = CON.m_continuation(DYN,S,AM,ST);
    end

end