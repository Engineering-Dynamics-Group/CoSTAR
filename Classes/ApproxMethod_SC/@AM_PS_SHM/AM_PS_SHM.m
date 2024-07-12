% Class Shoot provides functions which return the residual for shooting
% algorithms. Shoot is a subclass of ApproxMethod
classdef AM_PS_SHM < ApproxMethod

    properties
        odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);                    %options for the ODE integrators.     
        solver string = 'ode45';
        solver_function function_handle
        
        %Inherited Properties
        % res function_handle
    end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    methods(Static)                                                        %Static: Method can be called without creating an object of class SHM;

        s_PS_SHM_gatekeeper(GC,system,opt_approx_method,opt_init);          %Gatekeeper method, which is called by the static AM_gatekeeper method  
        help_text = s_help_opt_approx_method_PS_SHM();                        %Help text for the opt_approx_method option structure  
        help_text = s_help_opt_init_PS_SHM();                                 %Help text for the opt_init option structure  

    end
    %%%%%%%%%%%%%%%
    methods
        %% Constructor
        function obj = AM_PS_SHM(DYN)     
            obj = updateoptions(obj,DYN.opt_approx_method);                 %updateoptions method is a general method
            obj = obj.getIV(DYN);                                           %Set initial value (Has to be set here, because residual accesses iv
            obj = setSolver(obj,obj.solver);
        end
        
        %% Methods
        %interface methods
        obj = IF_up_res_data(obj,CON);                                     %This methods modifies the superclass method
        obj = getIV(obj,DYN);
   
        %% Functions for shooting algorithms
        f = SHM_single_fun(obj,y,DYN);                                   %Method required for SHM_single_nonauto
        f = SHM_single_auto_fun(obj,y,DYN);                              %Method required for SHM_single_auto


    end
end