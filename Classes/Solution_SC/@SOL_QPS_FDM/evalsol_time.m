% This function is a method of Solution subclass SOL_QPS_FDM and is called by the solget method of the superclass Solution
% The method evaluates time dependent solutions for the solution curve indices defined in options.index
% The calculated method solution vector s is interpolated by splines and is evaluated over the desired time intervall
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s:       Time solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @t:       Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function  [s,mu,t] = evalsol_time(obj,DYN,options)
    
    % Parameters
    dim = DYN.system.dim;                           % Dimension of state space
    index = options.index;                          % Indices of the solution curve where the solution is to be evaluated
    n_evals = length(index);                        % Number of evaluations
    res = options.resolution;                       % Desired resolution of output s (= number of discretized points in time over one period)

    if isfield(DYN.opt_approx_method,'n_int_1')
        n_int_1 = DYN.opt_approx_method.n_int_1;    % Get the number of intervals in theta_1-direction (n_int_2 can be calculated)
    else
        n_int_1 = 50;                               % Default value of n_int_1. Must be set again because the ApproxMethod object is not available here
    end
    if isfield(DYN.opt_approx_method,'n_int_2')
        n_int_2 = DYN.opt_approx_method.n_int_2;    % Get the number of intervals in theta_2-direction
    else
        n_int_2 = length(obj.s(:,index(1))) / (dim*n_int_1);    % Calculate the number of intervals in theta_2-direction
    end    


    % Create the time vector for the output. t(:,1,i) (i = 1,...,n_evals) defines the points in time at which the system states will be handed back after interpolation
    % In contrast to periodic solutions, the time vector is the same for each evaluation i. That is because there is no periodic time if the solution is quasi-periodic
    % Therefore, it does not make sense to calculate a specific time vector for each evaluation
    if isfield(options,'interval')                                              % If there is a user-defined time interval
        t(:,1) = linspace(options.interval(1), options.interval(2), res);       % User-defined time vector for output
    else
        t(:,1) = linspace(0, 2*pi, res);                                        % Default time vector for output (FGM and SHM: Upper limit is 2*pi*(1-1/res) -> can be changed?)
    end
    t = repmat(t,[1,1,n_evals]);                                                % Time vector is the same for each evaluation


    % Preallocate s and Z_i_interp
    s = zeros(res, dim, n_evals);                   % Preallocate s
    Z_i_interp = zeros(dim, n_int_1+1, n_int_2+1);  % Z_i_interp is preallocated because the array must already be present for the method of how it is "filled"


    % Calculate the theta coordinates which belong to the grid of Z_i_interp (needed for interpolation, see below)
    theta_1_interp = 0 : 2*pi/n_int_1 : 2*pi;       % Create the theta_1-coordinates ("x1"-values for the interpolation)
    theta_2_interp = 0 : 2*pi/n_int_2 : 2*pi;       % Create the theta_2-coordinates ("x2"-values for the interpolation)


    % Loop for the evaluations
    for i = 1:n_evals
        
        % The interpolation data Z_i_csape must be evaluated at t(:,1,i). However, the calculated state vectors in Z_i are defined on the 2D hypertime grid [0, 2*pi] x [0, 2*pi]
        % Thus, Z_i must be interpolated on the hypertime grid and also must be evaluated on the hypertime grid because the hypertime grid cannot be "rescaled" to a time grid
        % (In contrast to periodic solutions where the hypertime grid can easily be rescaled to a time grid because the scalar equation theta = omega * t holds)
        % In order to evaluate the interpolated data on the hypertime grid, hypertime coordinates, corresponding to the evaluation time points t(:,1,i), must be calculated
        % Due to the periodic conditions z(theta_1,theta_2) = z(theta_1+2*pi,theta_2) and z(theta_1,theta_2) = z(theta_1,theta_2+2*pi), the mentioned hypertime coordinates are equivalent to mod(hypertime_coordinates,2*pi)
        % This has the advantage that the state vectors only need to be interpolated in the [0, 2*pi] x [0, 2*pi] grid, since the values of mod(hypertime_coordinates,2*pi) are always on that grid
        Omega_i = obj.freq(:,index(i));             % Get the current frequencies (if system is partly autonomous: omega_2 is the autonomous frequency)
        theta_1_i_start = Omega_i(1)*t(1,1,i);      % theta_1-value at bottom limit of time interval
        theta_1_i_end = Omega_i(1)*t(end,1,i);      % theta_1-value at upper limit of time interval
        theta_2_i_start = Omega_i(2)*t(1,1,i);      % theta_2-value at bottom limit of time interval
        theta_2_i_end = Omega_i(2)*t(end,1,i);      % theta_2-value at upper limit of time interval
        theta_1_i_eval = mod(linspace(theta_1_i_start, theta_1_i_end, res), 2*pi);      % Create the hypertime coordinates theta_1_eval for the evaluation
        theta_2_i_eval = mod(linspace(theta_2_i_start, theta_2_i_end, res), 2*pi);      % Create the hypertime coordinates theta_2_eval for the evaluation
        
        % Get the current solution vector and reshape it to an array of dimension [dim x n_int_1 x n_int_2]
        s_i = obj.s(:,index(i));                    % Get the current solution vector s of the approximation method
        Z_i = reshape(s_i,[dim,n_int_1,n_int_2]);   % Reshape the solution vector (this structure is needed by csape)
        
        % The state vectors Z_i must be interpolated on the 2D hypertime grid [0, 2*pi] x [0, 2*pi]. However, Z_i covers the grid [0, 2*pi-DeltaTheta_1] x [0, 2*pi-DeltaTheta_2]
        % Therefore, a new matrix Z_i_interp covering the grid [0, 2*pi] x [0, 2*pi] is created ("y"-values for the interpolation)
        % The corresponding hypertime vectors were already created above since they are the same for each evaluation i
        Z_i_interp(:,1:end-1,1:end-1) = Z_i;        % Create a new array which will include the state vectors at (theta_1=2*pi, theta_2) and (theta_1, theta_2=2*pi)
        Z_i_interp(:,end,1:end-1) = Z_i(:,1,:);     % Add z(theta_1=2*pi, theta_2) = z(theta_1=0, theta_2) (for all theta_2 and for all dimensions of state space)
        Z_i_interp(:,1:end-1,end) = Z_i(:,:,1);     % Add z(theta_1, theta_2=2*pi) = z(theta_1, theta_2=0) (for all theta_1 and for all dimensions of state space)
        Z_i_interp(:,end,end) = Z_i(:,1,1);         % Add z(theta_1=2*pi, theta_2=2*pi) = z(theta_1=0, theta_2=0) (for all dimensions of state space)

        % Interpolate the state vectors on the 2D hypertime grid using csape
        % Here, csape uses periodic end conditions, i.e. the first and second derivatives at the left end are matched with those at the right end (for each grid dimension)
        % To interpolate vectored grid data, the vectors must be stored in a 3-dimensional array with Z(:,i,j) being the vector at the grid indices (i,j)
        Z_i_csape = csape({theta_1_interp,theta_2_interp},Z_i_interp,{'periodic','periodic'});
        Z_i_eval = fnval(Z_i_csape,{theta_1_i_eval,theta_2_i_eval});        % Evaluate the interpolation at the evaluation coordinates {theta_1_i_eval,theta_2_i_eval}

        % fnval has evaluated the interpolated hypertime map for all theta_1_i_eval and theta_2_i_eval, i.e. the output Z_i_eval is a map again ( size(Z_i_eval) = [dim x res x res] )
        % For time evaluation, only the state vectors on the main diagonal of Z_i_map(j,:,:) (j = 1,...,dim) are needed
        s_i_eval = Z_i_eval(:,1:res+1:end);             % This extracts the state vectors from the main diagonal of Z_i_eval(j,:,:)
        s(:,:,i) = s_i_eval';                           % Transponse s_i_eval since size(s_i_eval) = [dim x res]

    end


    % Get the mu values
    mu = obj.mu(index);                         % It is ensured by S.solget_up_index that the elements of options.index are unambigious 

end