% This function is a method of Solution subclass SOL_PS_FDM
% It is called by the solget method of the superclass Solution
% The method evaluates time dependent solutions for the solution curve indices defined in options.index
% The calculated method solution vector s is interpolated by splines and is evaluated in the desired time intervall
% ATTENTION: This method must be adapted when error control for finite differences is implemented
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
    dim = DYN.system.dim;                       % Dimension of state space
    index = options.index;                      % Indices of the solution curve where the solution is to be evaluated
    n_evals = length(index);                    % Number of evaluations
    res = options.resolution;                   % Desired resolution of output s (= number of discretized points in time over one period)


    % Preallocate t and s (t is calculated in the loop as it can change for every evaluation due to periodic time T being different)
    t = zeros(res, 1, n_evals);
    s = zeros(res, dim, n_evals);


    % Loop for the evaluations
    for i = 1:n_evals

        T_i = 2*pi/obj.freq(1,index(i));            % Get the current periodic time
        
        % Create the time vector for the output. t(:,1,i) defines the points in time at which the system states will be handed back after interpolation
        if isfield(options,'interval')                                              % If there is a user-defined time interval 
            t(:,1,i) = linspace(options.interval(1), options.interval(2), res);     % User-defined time vector for output
        else
            t(:,1,i) = linspace(0, T_i, res);                                       % Default time vector for output
        end
        
        % The interpolation data Z_i_csape must be evaluated at t(:,1,i)
        % Due to periodic condition z(t) = z(t+T), an evaluation at t_i_eval = mod(t(:,1,i),T) is equivalent to an evaluation at t(:,1,i)
        % This has the advantage that the state vectors only need to be interpolated in the interval [0, T], since the values of t_i_eval are always in the interval [0, T]
        t_i_eval = mod(t(:,1,i),T_i);

        % Get the current solution vector s and reshape it to a matrix Z_i
        s_i = obj.s(:,index(i));                    % Get the current solution vector s_i of the approximation method
        Z_i = reshape(s_i,dim,[]);                  % Reshape s_i to Z_i = [z(0), ... , z(T-DeltaT)] -> size(Z_i) = [dim n_int] 
        
        % The state vectors in Z_i must be interpolated in the time interval [0, T]. However, Z_i covers the interval [0, T-DeltaT]
        % Therefore, a new matrix Z_i_interp = [z(0), ..., z(T)] = [Z_i, Z_i(:,1)] covering the time interval [0, T] is created ("y"-values for the interpolation)
        % Additionally, a corresponding time vector needs to be created ("x"-values for the interpolation)
        Z_i_interp = [Z_i, Z_i(:,1)];                                   % Add z(t=T) = z(t=0) (periodic condition) to cover the time interval [0, T]
        t_i_interp = linspace(0, T_i, size(Z_i_interp,2));              % Create the time vector belonging to Z_i_interp (needed for interpolation)

        % Interpolation
        % s(:,:,i) = interp1(t_i_interp, Z_i_interp', t_i_eval, 'spline');   % Interpolate the state vectors in each dimension
        % The structure of Z_i_interp' is needed because interp1 interpolates each column of Z_i_interp' (= each system state) seperately (see also MATLAB doc.)
        % This has the positive side effect that the output of interp1 already has the correct size for s(:,:,i)
        % size(t_i_interp)  = [1         (n_int+1)] 
        % size(Z_i_interp') = [(n_int+1)    dim   ]     -->   size(s(:,:,i)) = [res  dim]
        % size(t_i_eval)    = [1            res   ] 
        % NEW:
        % Using csape instead of interp1. csape interpolates the state vectors in each dimension using periodic end conditions
        % This matches the first and second derivatives at the left end with those at the right end)
        % In contrast to interp1, csape needs the "y"-data to be an [dim x (n_int+1)] array
        Z_i_csape = csape(t_i_interp, Z_i_interp, 'periodic');          % Interpolation. Output of csape is a struct (piecewise polynomial form)   
        s(:,:,i) = fnval(Z_i_csape, t_i_eval)';                         % Evaluate the interpolated curve at t(:,1,i). The dimension of the output is [dim x res]

    end


    % Get the mu values
    mu = obj.mu(index);                         % It is ensured by S.solget_up_index that the elements of options.index are unambigious 

end