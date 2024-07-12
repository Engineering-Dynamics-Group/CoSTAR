%Class DynamicalSystem creates an object which contains all relevant
%informations about the dynamical system. The ODE must be given at the
%moment as either an algebraic equation (equilibrium) 
% 
% 0 = f(z,p)
%
% or a first order ODE
%
%z' = f(t,z,p)
%
%where p is a column-vector of parameter's of the system 

classdef DynamicalSystem 
    
    properties
        info string                                                        %information or/and description of the system    

        %system structure
        order   uint8                                                      %order of the ode
        rhs     function_handle                                            %right hand side of the ode
        dim     double                                                     %dimension of the state space
        param = {0};                                                       %array for the parameter vector
        
        %opt_sol structure
        sol_type        string                                             %equilibrium, periodic, quasi-periodic
        approx_method   string                                             %shooting, fd, fgm
        cont            string                                             %indicates if a continuation should be done or not
        stability       string                                             %indicates if stability should be computed or not

        act_param = 1;                                                     %Defines the index of the continuation parameter in the param array
        non_auto_freq = []                                                 %contains all non-autonomous frequencies (this may be type function handle)
        auto_freq = []                                                     %contains all initial values for the autonomous frequencies

           %Internally set
        n_auto uint32                                                      %Gives the number of autonomous frequencies
        n_freq                                                             %Number of Base Frequencies
        DYN_id  string                                                     %identification string, which is used to match DynamicalSystem object and Solution object
    end

    properties(SetAccess = private)
        %% All options structures are stored here for evaluation purposes
        system struct                                                      %struct defining the system
        opt_sol struct                                                     %options (struct) for the solution branch
        opt_approx_method struct                                           %options (struct) for approximation method
        opt_init struct
        opt_cont struct                                                    %options (struct) for continuation
        opt_stability struct                                               %options (struct) for stability computation
    end                                                %identification string, which is used to match DynamicalSystem object and Solution object


    
    methods
        %% Constructor
        function obj = DynamicalSystem(options)                            %Constructor
            obj = updateoptions(obj,options.system);                       %Assign all variables to  the properties   
            obj = updateoptions(obj,options.opt_sol);                      %Assign all variables to  the properties   

            % Options are saved here for documentation purposes
            obj.system              = options.system;
            obj.opt_sol             = options.opt_sol;
            obj.opt_approx_method   = options.opt_approx_method;
            obj.opt_init            = options.opt_init;
            obj.opt_cont            = options.opt_cont;
            obj.opt_stability       = options.opt_stability;
            
            % If user supplied 'eq', 'ps', 'qps', 'fgm', 'fdm' or 'shm': change to "long" names since code is written using them
            if strcmpi(obj.sol_type,'eq');      obj.sol_type = 'equilibrium';
            elseif strcmpi(obj.sol_type,'ps');  obj.sol_type = 'periodic';
            elseif strcmpi(obj.sol_type,'qps'); obj.sol_type = 'quasiperiodic';
            end
            if strcmpi(obj.approx_method,'fgm');        obj.approx_method = 'fourier-galerkin';
            elseif strcmpi(obj.approx_method,'shm');    obj.approx_method = 'shooting';
            elseif strcmpi(obj.approx_method,'fdm');    obj.approx_method = 'finite-difference';
            end

            obj = obj.updatefreq();
            
            obj.DYN_id = strjoin(['DYN','ID',cellstr(string(clock))],'_'); %Creates a unique ID: This ID is used to match Dynamical systems and solution objects in postprocessing.
        end
        
        %% Methods
        obj = updatefreq(obj,opt_init);                                             %This function sets the frequency and intial conditions vectors correctly.
    end

    methods(Static)

        s_DYN_gatekeeper(GC,system,opt_sol)                                  %Gatekeeper for the DynamicalSystem Class
        help_struct = s_help_system();                                                     %help file for the system costaropts structure    
        help_struct = s_help_opt_sol();                                                    %help file for the opt_sol costaropts structure    


    end
    
end




