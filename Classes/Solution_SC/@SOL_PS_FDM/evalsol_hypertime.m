% This function is a method of Solution subclass SOL_PS_FDM and is called by the solget method of the superclass Solution
% It calculates the hypertime manifolds of a solution, which is 1-dimensional
% The method evaluates hypertime solutions for the solution curve indices defined in options.index
% ATTENTION: This method must be adapted when error control for finite differences is implemented
%
% @obj:       Solution subclass object
% @DYN:       DynamicalSystem object
% @options:   options structure for postprocessing solutions
%
% @s:         Hypertime solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:        Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @hypertime: Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:    Number of curve points to be evaluated 

function [s,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

    % Parameters
    dim = DYN.system.dim;                           % Dimension of state space
    index = options.index;                          % Indices of the solution curve where the solution is to be evaluated
    n_evals = length(index);                        % Number of evaluations
    res = options.resolution;                       % Desired resolution of output s (= number of discretized points in time over one period)

    if isfield(DYN.opt_approx_method,'n_int')
        n_int = DYN.opt_approx_method.n_int;        % Get the number of intervals in theta-direction
    else
        n_int = length(obj.s(:,index(1))) / dim;    % Calculate the number of intervals if user did not supply n_int (the ApproxMethod object is not available here)
    end


    % Calculate hypertime array for the output
    hypertime(:,1) = linspace(0, 2*pi, res);        % Create the first two dimensions of the hypertime array. First array dimension defines the hypertime coordinates of the output
    hypertime = repmat(hypertime,[1,1,n_evals]);    % The hypertime coordinates are the same for each evaluation
    

    % Preallocate s
    s = zeros(res, dim, n_evals);


    % Create the theta coordinates which belong to the state vectors stored in Z_i_interp (needed for interpolation, see below)
    theta_interp = 0 : 2*pi/n_int : 2*pi;           % Create the theta coordinates ("x"-values for the interpolation)


    % Loop for the evaluations
    for i = 1:n_evals

        % Get the current solution vector and reshape it to a matrix
        s_i = obj.s(:,index(i));                    % Get the current solution vector s_i of the approximation method
        Z_i = reshape(s_i,dim,[]);                  % Reshape s_i to Z_i = [z(0), ... , z(2*pi-DeltaTheta)] -> size(Z_i) = [dim n_int]

        % The state vectors in Z_i must be interpolated in the theta interval [0, 2*pi]. However, Z_i covers the interval [0, 2*pi-DeltaTheta]
        % Therefore, a new matrix Z_i_interp = [z(0), ..., z(2*pi)] = [Z_i, Z_i(:,1)] covering the theta interval [0, 2*pi] is created ("y"-values for the interpolation)
        % The corresponding hypertime vector was already created above since it is the same for each evaluation i
        Z_i_interp = [Z_i, Z_i(:,1)];               % Add z(theta=2*pi) = z(theta=0) (periodic condition) to cover the output interval [0, 2*pi]

        % Interpolate the state vectors in each dimension using csape
        % Here, csape uses periodic end conditions, i.e. the first and second derivatives at the left end are matched with those at the right end
        Z_i_csape = csape(theta_interp, Z_i_interp, 'periodic');       % Interpolation. Output of csape is a struct (piecewise polynomial form)
        s(:,:,i) = fnval(Z_i_csape, hypertime(:,1,i))';                % Evaluate the interpolated curve at the evaluation coordinates. The dimension of the output is [dim x res]
    
    end


    % Get the mu values
    mu = obj.mu(index);                             % It is ensured by S.solget_up_index that the elements of options.index are unambigious 


    % Old code: evalsol_time is called and the time array is rescaled to [0, 2*pi] afterwards
    % This method, however, needs a bit more computation than the code above
    %{
    [s,mu,time] = obj.evalsol_time(DYN,options);                        % Get the time dependent solutions

    % Normalize time intervall [0, 2*pi/freq] to hypertime intervall [0, 2*pi]: Each time vector t(:,1,i) is multiplied by obj.freq(1,options.index(i))
    % hypertime = zeros(size(t));
    % for i = 1:length(options.index)
    %     hypertime(:,1,i) = time(:,1,i) .* obj.freq(1,options.index(i));
    % end
    % Or, in short:
    hypertime = time.*permute(obj.freq(1,options.index),[1,3,2]);       % size(permute(obj.freq(1,options.index),[1,3,2])) = [1  1  length(options.index)]
    %}

end