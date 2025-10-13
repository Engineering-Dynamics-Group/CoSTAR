% Class DynamicalSystem creates an object which contains all relevant information about the dynamical system
% The ODE must be given either as an algebraic equation (equilibrium) 0 = f(z,param) or as a first order ODE
% z' = f(t,z,param), where param is a cell array of parameters of the system

classdef DynamicalSystem 
    
    %% Properties    
    properties(Constant)
        % Constant properties cannot be changed (read-only)
        costar_version = '3.1.3';                                       %release versions (e.g. 3.1) should be defined as ' 3.1 ' (see the space between the ' and the numbers)
    end


    properties
        info string                                                     %information or/and description of the system    

        %system structure
        order   uint8                                                   %order of the ode
        rhs     function_handle                                         %right hand side of the ode
        dim     double                                                  %dimension of the state space
        param = {0};                                                    %array for the parameter vector
        
        %opt_sol structure
        sol_type        string                                          %equilibrium, periodic, quasi-periodic
        approx_method   string                                          %shooting, fd, fgm
        cont            string                                          %indicates if a continuation should be done or not
        stability       string                                          %indicates if stability should be computed or not

        act_param = 1;                                                  %Defines the index of the continuation parameter in the param array
        non_auto_freq = []                                              %contains all non-autonomous frequencies (this may be type function handle)
        auto_freq = []                                                  %contains all initial values for the autonomous frequencies
        freq_limit = 1e-4;                                              %CoSTAR stops when the base frequency(s) fall below this limit (ATTENTION: Also set in s_DYN_gatekeeper)

        %Internally set
        n_auto uint32                                                   %Gives the number of autonomous frequencies
        n_freq                                                          %Number of Base Frequencies
        DYN_id  string                                                  %identification string, which is used to match DynamicalSystem object and Solution object
    end


    properties(SetAccess = private)
        % All options structures are stored here for evaluation purposes
        system struct                                                   %struct defining the system
        opt_sol struct                                                  %options (struct) for the solution branch
        opt_approx_method struct                                        %options (struct) for approximation method
        opt_init struct
        opt_cont struct                                                 %options (struct) for continuation
        opt_stability struct                                            %options (struct) for stability computation
    end


    properties(Hidden)
        % Hidden properties are not visible in property lists or in results from calls to get, set, or the properties functions
        create_log = true;                                              %defines whether a log file is created (used by the development team)
        display = 'iter';                                               %Controls the command window output
    end



    %% Methods
    methods

        % Constructor
        function obj = DynamicalSystem(options)                         %Constructor
            obj = updateoptions(obj,options.system);                    %Assign all variables to  the properties   
            obj = updateoptions(obj,options.opt_sol);                   %Assign all variables to  the properties   

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
            
            % Creates an unique ID (used to match DynamicalSystem and Solution objects in postprocessing)
            obj.DYN_id = append('ID_', char(datetime('now','TimeZone','local','Format','yyyy-MM-dd_HH.mm.ss')));  

            lastwarn('');                                               %Reset lastwarn (needed for the "Finished" message)        
        end
        

        % Other methods
        obj = updatefreq(obj,opt_init);                                 %This function sets the frequency and initial conditions vectors correctly.
    
    end


    % Static methods
    methods(Static)

        s_DYN_gatekeeper(GC,system,opt_sol)                             %Gatekeeper for the DynamicalSystem Class
        help_struct = s_help_system();                                  %help file for the system costaropts structure    
        help_struct = s_help_opt_sol();                                 %help file for the opt_sol costaropts structure    

    end
    

end