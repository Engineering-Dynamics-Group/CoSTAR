% Class AM_PS_FDM provides functions which return the residual for finite-difference method algorithms (perdiodic solution)
% AM_PS_FDM is a subclass of ApproxMethod

classdef AM_PS_FDM < ApproxMethod

    properties                  % Define the properties, which are supplied in the opt_approx_method and opt_init sturctures
        
        % ATTENTION: If the default values are changed, they must be changed in s_PS_FDM_gatekeeper as well!
        n_int = 100;            % Number of hyper-time intervals DeltaTheta into which the hyper-time period 2*pi is divided, i.e. 2*pi = n_int * DeltaTheta
        scheme = 'central';     % Discretization scheme used to approximate the derivation dz(theta_i)/dtheta
        approx_order = 6;       % Order of finite-difference approximation of dz(theta_i)/dtheta
        points                  % Stores the local grid point indices sigma_k which are used to approximate dz/dtheta at theta_i using z_(i+sigma_k) = z(theta_i + sigma_k * DeltaTheta)               
        c0                      % Fourier-coefficient of 0-th order to create the initial solution vector
        c1                      % Fourier-coefficient of 1-st order cosine term to create the initial solution vector
        s1                      % Fourier-coefficient of 1-st order sine term to create the initial solution vector
        fdm_sol                 % Already calculated method solution vector s used as an initial value to calculate the initial solution

        % Inherited properties from superclass ApproxMethod
        % res function_handle                                           % Residual function defined by a solution method (e.g. Shoot)
        % iv
        % error_control = 'off';                                        % Switch for using the error_control feature (default value: off)
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = private)
    
        p_weights;              % Stores the weighting factors needed to approximate dz(theta_i)/dtheta by dz_i/dtheta = 1/DeltaTheta * sum_(k=1)^p ( w_(sigma_k) * z_(i+sigma_k) )
        p_w_mat_J;              % Stores a diagonal sparse matrix containing the weights according to a certain pattern. Needed for the Jacobian matrix
        p_ind_blkdiag_mat;      % Stores the indices which are needed to create a block diagonal sparse matrix. This accelerates the calculation of the Jacobian matrix

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static)                                                     % Static methods: can be called without creating an object of class AM_PS_FDM

        s_PS_FDM_gatekeeper(GC,system,opt_approx_method,opt_init);      % Gatekeeper method, which is called by the static s_AM_gatekeeper method to check the inputs at the beginning 
        
        help_struct = s_help_opt_approx_method_PS_FDM();                % Help struct for the opt_approx_method option structure 
        help_struct = s_help_opt_init_PS_FDM();                         % Help struct for the opt_init option structure

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        %% Constructor
        % @DYN: DynamicalSystem class object

        function obj = AM_PS_FDM(DYN)     
            
            obj = updateoptions(obj,DYN.opt_approx_method);             % Updates the properties with the values set by the user in struct opt_approx_method
            obj = updateoptions(obj,DYN.opt_init);                      % Updates the properties with the values set by the user in opt_init
            
            obj = getWeights(obj,DYN);                                  % Get (or calculate) the weights used to approximate dz(theta_i)/dtheta

            obj = getIV(obj,DYN);                                       % Set initial value (must be set here, because function initial_solution needs it)
        
        end
        
        %% Methods
        % @DYN: DynamicalSystem class object
        % @CON: Continuation class object
        % @y:   solution vector of nonlinear equation system
        % @y0:  initial value of solution vector of nonlinear equation system
        
        obj = getWeights(obj,DYN);              % Method to get (or calculate) the weights used to approximate dz(theta_i)/dtheta
        obj = getIV(obj,DYN);                   % Method that generates an initial value for nonlinear equation solver (e.g. fsolve) to start from. Uses options defined in opt_init structure

        [F,J] = corr_fun_init_FDM(obj,y,y0);        % Method that provides the residuum vector function and the corresponding Jacobian matrix for fsolve to calculate the initial solution
        [F,J] = corr_fun_FDM(obj,y,CON);            % Method that provides the residuum vector function and the corresponding Jacobian matrix for fsolve during continuation
        [res,J_res] = PS_FDM_residuum(obj,y,DYN);   % Method that builds the residuum of the finite-difference equation system and its Jacobian matrix
        
        obj = IF_up_res_data(obj,CON);          % Interface method: Used to pass information between continuation algorithm an this subclass
        IC  = getIC(obj,y,DYN,n_shoot);         % Method that extracts the state space vector z(theta=0), which is needed for stability calculation via the shooting method

    end

end