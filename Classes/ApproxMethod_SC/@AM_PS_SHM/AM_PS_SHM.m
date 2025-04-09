% Class Shoot provides functions which return the residual for shooting
% algorithms. Shoot is a subclass of ApproxMethod
classdef AM_PS_SHM < ApproxMethod

    properties
        odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);                 % options for the ODE integrator     
        solver string = 'ode45';
        solver_function function_handle
        n_shoot = 2;                                                    % Number of shooting points for multiple shooting 

        y_old                                                           % "Dummy" property so that the if condition in m_continuation, line 78, can be used for PS and QPS
        
        %Inherited Properties
        % res function_handle
    end

    %%%%%%%%%%%%%%%

    methods(Static)                                                     % Static: Method can be called without creating an object of class SHM;

        s_PS_SHM_gatekeeper(GC,system,opt_approx_method,opt_init);      % Gatekeeper method, which is called by the static AM_gatekeeper method  
        help_text = s_help_opt_approx_method_PS_SHM();                  % Help text for the opt_approx_method option structure  
        help_text = s_help_opt_init_PS_SHM();                           % Help text for the opt_init option structure  

    end

    %%%%%%%%%%%%%%%

    methods

        % Constructor
        function obj = AM_PS_SHM(DYN)     
            obj = updateoptions(obj,DYN.opt_approx_method);             % updateoptions method is a general method
            obj = setSolver(obj,obj.solver);
            obj = obj.getIV(DYN);                                       % Set initial value (Has to be set here, because residual accesses iv
        end
        
        % Interface Methods
        obj = IF_up_res_data(obj,CON,DYN);                              % This methods modifies the superclass method
        obj = getIV(obj,DYN);
   
        % Methods for shooting algorithms
        [res,J_res] = PS_SHM_residuum(obj,y,DYN);                       % Residuum function
        [F,J] = fun_Jac_wrapper(obj,y,CONT);                            % Function wrapper for fsolve to evaluate jacobian "analytically"
        [F,J] = fun_Jac_wrapper_init(obj,y,y0);                         % Function wrapper_init for fsolve to evaluate jacobian "analytically" for initial solution

    end
    
end