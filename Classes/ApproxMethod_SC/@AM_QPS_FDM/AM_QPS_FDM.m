% Class AM_QPS_FDM provides functions which return the residual for finite-difference method algorithms (quasi-perdiodic solution)
% AM_QPS_FDM is a subclass of ApproxMethod

classdef AM_QPS_FDM < ApproxMethod

    properties                  % Define the properties, which are supplied in the opt_approx_method and opt_init sturctures
        
        % ATTENTION: If the default values are changed, they must be changed in s_QPS_FDM_gatekeeper as well!
        n_int_1 = 50;           % Number of hyper-time intervals DeltaTheta_1 into which the hyper-time period 2*pi is divided in theta_1-direction, i.e. 2*pi = n_int_1 * DeltaTheta_1 
                                % -> n_int_1 ALSO NEEDS TO BE SET IN evalsol_hypertime and evalsol_time (because the AM-object is not available there) !
        n_int_2 = 50;           % Number of hyper-time intervals DeltaTheta_2 into which the hyper-time period 2*pi is divided in theta_2-direction, i.e. 2*pi = n_int_2 * DeltaTheta_2
        scheme_1 = 'central';   % Discretization scheme used to approximate the derivation dz(theta_1_i,theta_2_j)/dtheta_1
        scheme_2 = 'central';   % Discretization scheme used to approximate the derivation dz(theta_1_i,theta_2_j)/dtheta_2
        approx_order_1 = 6;     % Order of finite-difference approximation of dz(theta_1_i,theta_2_j)/dtheta_1
        approx_order_2 = 6;     % Order of finite-difference approximation of dz(theta_1_i,theta_2_j)/dtheta_2
        points_1;               % Stores the local grid point indices sigma_1_k which are used to approximate dz(theta_1_i,theta_2_j)/dtheta_1 at using z_(i+sigma_1_k)_j = z(theta_i + sigma_1_k * DeltaTheta_1, theta_j)
        points_2;               % Stores the local grid point indices sigma_2_k which are used to approximate dz(theta_1_i,theta_2_j)/dtheta_2 at using z_i_(j+sigma_2_k) = z(theta_i, theta_j + sigma_2_k * DeltaTheta_2)
        weights_1;              % Stores the weighting factors needed to approximate dz(theta_1_i,theta_2_j)/dtheta_1 by dz_i_j/dtheta_1 = 1/DeltaTheta_1 * sum_(k=1)^p ( w_1_(sigma_1_k) * z_(i+sigma_1_k)_j )
        weights_2;              % Stores the weighting factors needed to approximate dz(theta_1_i,theta_2_j)/dtheta_2 by dz_i_j/dtheta_2 = 1/DeltaTheta_2 * sum_(k=1)^p ( w_2_(sigma_2_k) * z_i_(j+sigma_2_k) )
        c0                      % Fourier-coefficient of 0-th order to create the initial solution vector
        c1_matrix               % Fourier-coefficients (stored in a matrix) of 1-st order cosine terms to create the initial solution vector
        s1_matrix               % Fourier-coefficients (stored in a matrix) of 1-st order sine terms to create the initial solution vector
        fdm_sol                 % Already calculated method solution vector s used as an initial value to calculate the initial solution
        n_int_1_fdm_sol         % Number of intervals n_int_1 corresponding to fdm_sol. Needs to be provided by user if it is different from n_int_1
        n_int_2_fdm_sol         % Number of intervals n_int_2 corresponding to fdm_sol. Needs to be provided by user if it is different from n_int_2

        % Inherited properties from superclass ApproxMethod
        % res function_handle                                           % Residual function defined by a solution method (e.g. Shoot)
        % iv
        % error_control = 'off';                                        % Switch for using the error_control feature (default value: off)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%

    properties(Access = private)
        p_w_1_mat_J;            % Stores a diagonal sparse matrix containing the weights_1 according to a certain pattern. Needed for the Jacobian matrix
        p_w_2_mat_J;            % Stores a diagonal sparse matrix containing the weights_2 according to a certain pattern. Needed for the Jacobian matrix
        p_ind_blkdiag_mat;      % Stores the indices which are needed to create a block diagonal sparse matrix. This accelerates the calculation of the Jacobian matrix

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static)                                                     % Static methods: can be called without creating an object of class AM_QPS_FDM

        s_QPS_FDM_gatekeeper(GC,system,opt_approx_method,opt_init);     % Gatekeeper method, which is called by the static s_AM_gatekeeper method to check the inputs at the beginning 
        
        help_struct = s_help_opt_approx_method_QPS_FDM();               % Help struct for the opt_approx_method option structure 
        help_struct = s_help_opt_init_QPS_FDM();                        % Help struct for the opt_init option structure

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        %% Constructor
        % @DYN: DynamicalSystem class object

        function obj = AM_QPS_FDM(DYN)     
            
            obj = updateoptions(obj,DYN.opt_approx_method);             % Updates the properties with the values set by the user in struct opt_approx_method
            obj = updateoptions(obj,DYN.opt_init);                      % Updates the properties with the values set by the user in opt_init
            
            obj = getWeights(obj,DYN);                                  % Get (or calculate) the weights used to approximate dz_i_j/dtheta

            obj = getIV(obj,DYN);                                       % Set initial value (must be set here, because function initial_solution needs it)

        end
        
        %% Methods
        % @DYN: DynamicalSystem class object
        % @CON: Continuation class object
        % @y:   solution vector of nonlinear equation system
        % @y0:  initial value of solution vector of nonlinear equation system
        
        obj = getWeights(obj,DYN);              % Method to get (or calculate) the weights used to approximate dz_i_j/dtheta
        obj = getIV(obj,DYN);                   % Method that generates an initial value for nonlinear equation solver (e.g. fsolve) to start from. Uses options defined in opt_init structure

        [F,J] = corr_fun_init_FDM(obj,y,y0);        % Method that provides the residuum vector function and the corresponding Jacobian matrix for fsolve to calculate the initial solution
        [F,J] = corr_fun_FDM(obj,y,CON);            % Method that provides the residuum vector function and the corresponding Jacobian matrix for fsolve during continuation
        [res,J_res] = QPS_FDM_residuum(obj,y,DYN);  % Method that builds the residuum of the finite-difference equation system and its Jacobian matrix
        
        obj = IF_up_res_data(obj,CON);          % Interface method: Used to pass information between continuation algorithm an this subclass
        % IC  = getIC(obj,y,DYN);               % Method that extracts the state space vector z(0,0) (periodic solution: needed needed for stability calculation via the shooting method)

    end

end