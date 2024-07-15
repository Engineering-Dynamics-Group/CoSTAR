% Class AM_QPS_SHM provides functions which return the residual for
% quasi-periodic shooting algorithms. SHM is a subclass of ApproxMethod
classdef AM_QPS_SHM < ApproxMethod

    properties
        ic                                                                  %Initial condition
        mu0                                                                 %Initial value of mu
        n                                                                   %Dimension of state-space
        n_char = 100;                                                       %Number of characteristics
        Ik                                                                  %Integration interval
        phi
        odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);                     %options for the ODE integrators.
        solver string = 'ode45';                                            %Ode-solver
        solver_function function_handle                                     %Ode-solver function
        
        tinit = 10000;
        deltat = 15000;
        dt = 0.1;

        y_old                                                               %Initial conditions of last solution point (for phase condition in autonomous case)
        Y_old                                                               %Values of manifold of last solution point (for phase condition in autonomous case)
        reso_phase = 50;                                                    %Resolution for time integration to determine phase condition
        
%         phase_shift
        
        
        %Inherited Properties
        % res function_handle
    end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    methods(Static)                                                                 % Static: Method can be called without creating an object of class SHM;

        s_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init);       % Gatekeeper method, which is called by the static ST_gatekeeper method
        help_text = s_help_opt_approx_method_QPS_SHM();                             % Help text for the opt_approx_method option structure
        help_text = s_help_opt_init_QPS_SHM();                                      % Help text for the opt_init option structure

    end
    %%%%%%%%%%%%%%%
    methods
        %% Constructor
        function obj = AM_QPS_SHM(DYN)
            obj = updateoptions(obj,DYN.opt_approx_method);                         % updateoptions method is inherited from SolutionType
            obj = setSolver(obj,obj.solver);                                        % set ode Solver for time integration
            obj.n = DYN.dim;                                                        % Get dimension of system
            obj.phi = zeros(obj.n,1);                                               % Set phase for characteristics to zero
            
            if(isfield(DYN.opt_init,'tinit')); obj.tinit = DYN.opt_init.tinit; end
            if(isfield(DYN.opt_init,'deltat')); obj.deltat = DYN.opt_init.deltat; end
            if(isfield(DYN.opt_init,'dt')); obj.dt = DYN.opt_init.dt; end
            
            obj = getIV(obj,DYN);                                                   % Get initial value for starting solution
        end

        %% Methods
        %interface methods
        obj = IF_up_res_data(obj,CON,DYN,S);                                %This methods modifies the superclass method
        obj = getIV(obj,DYN);                                               %Method generates initial value for shooting method from point in state-space

        %% Methods for shooting algorithms
        [f,J] = qp_SHM_non_auto_fun(obj,y,DYN);                           % Residual for non-autonomous case for quasi-periodic shooting algorithm
        [f,J] = qp_SHM_mixed_fun(obj,y,DYN);                              % Residual for mixed case for quasi-periodic shooting algorithm
        [f,J] = qp_SHM_auto_fun(obj,y,DYN);                               % Residual for full-autonomous case for quasi-periodic shooting algorithm

        [f,J] = fun_Jac_wrapper(obj,y,DYN,CONT);                            %Function wrapper for fsolve to evaluate jacobian "analytically"
        [f,J] = fun_Jac_wrapper_init(obj,y,y0,DYN);                         %Function wrapper_init for fsolve to evaluate jacobian "analytically" for initial solution

        f = FcnWrapperODE2(obj,t,z,Fcn,PHI);                                %Function wrapper for ode-solver                             
        f = FcnWrapperODE5(obj,t,z,Fcn,PHI);                                %Function wrapper for ode-solver  
        
        P = poincare_int(obj,F,F1,Omega,Ik);                                %Poincare phase condition
    end
end