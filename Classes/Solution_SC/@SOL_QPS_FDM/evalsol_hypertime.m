% This function is a method of Solution subclass SOL_QPS_FDM and is called by the solget method of the superclass Solution
% It calculates the hypertime manifolds of a solution, which are 2-dimensional
% The method evaluates hypertime solutions for the solution curve indices defined in options.index
% ATTENTION: This method must be adapted when error control for finite differences is implemented
%
% @obj:       Solution subclass object
% @DYN:       DynamicalSystem object
% @options:   options structure for postprocessing solutions
%
% @s:         Hypertime solution array: This must(!) be a [options.resolution(1) x options.resolution(2) x state_space_dimension x n_evals] dimensional array !!!
% @mu:        Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @hypertime: Array of the time points: This must(!) be a [options.resolution(1) x options.resolution(2) x 2 x n_evals] dimensional array !!!
% n_evals:    Number of curve points to be evaluated 

function [s,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

    % Parameters
    dim = DYN.system.dim;                           % Dimension of state space
    index = options.index;                          % Indices of the solution curve where the solution is to be evaluated
    n_evals = length(index);                        % Number of evaluations
    if size(options.resolution,2) == 2
        res_1 = options.resolution(1);              % Desired resolution of output s in theta_1-direction if user defined [1x2] array
        res_2 = options.resolution(2);              % Desired resolution of output s in theta_2-direction if user defined [1x2] array
    else
        res_1 = options.resolution;                 % Desired resolution of output s in theta_1-direction if user defined scalar
        res_2 = res_1;                              % Desired resolution of output s in theta_2-direction if user defined scalar
    end
    
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


    % Calculate the hypertime array for the output
    theta_eval_1 = linspace(0, 2*pi, res_1);        % Coordinates of the output in theta_1-direction
    theta_eval_2 = linspace(0, 2*pi, res_2);        % Coordinates of the output in theta_2-direction
    [Theta_2, Theta_1] = meshgrid(theta_eval_2,theta_eval_1);     % Theta_2 is the first output of meshgrid so that hypertime(:,:,2) has the desired structure
    hypertime(:,:,1) = Theta_1;                     % Theta_1 values as a matrix. First row is the zero vector, second row only consists of values 2*pi/res etc.
    hypertime(:,:,2) = Theta_2;                     % Theta_2 values as a matrix. First column is the zero vector, second column only consists of values 2*pi/res etc.
    hypertime = repmat(hypertime,[1,1,1,n_evals]);  % The matrices are the same for each evaluation


    % Preallocate s and Z_i_interp
    s = zeros(res_1, res_2, dim, n_evals);          % Preallocate s
    Z_i_interp = zeros(dim, n_int_1+1, n_int_2+1);  % Z_i_interp is preallocated because the array must already be present for the method of how it is "filled"

    
    % Calculate the theta coordinates which belong to the grid of Z_i_interp (needed for interpolation, see below)
    theta_1_interp = 0 : 2*pi/n_int_1 : 2*pi;       % Create the theta_1-coordinates ("x1"-values for the interpolation)
    theta_2_interp = 0 : 2*pi/n_int_2 : 2*pi;       % Create the theta_2-coordinates ("x2"-values for the interpolation)


    % Loop for the evaluations
    for i = 1:n_evals

        % Get the current solution vector and reshape it to an array of dimension [dim x n_int_1 x n_int_2]
        s_i = obj.s(:,index(i));                    % Get the current solution vector s of the approximation method
        Z_i = reshape(s_i,[dim,n_int_1,n_int_2]);   % Reshape the solution vector (this structure is needed by csape)

        % The state vectors in Z_i must be interpolated on the 2D hypertime grid [0, 2*pi] x [0, 2*pi]
        % However, Z_i covers the grid [0, 2*pi-DeltaTheta_1] x [0, 2*pi-DeltaTheta_2]
        % Therefore, a new matrix Z_i_interp covering the grid [0, 2*pi] x [0, 2*pi] is created ("y"-values for the interpolation)
        % The corresponding hypertime vectors were already created above since they are the same for each evaluation i
        Z_i_interp(:,1:end-1,1:end-1) = Z_i;        % Create a new array which will include the state vectors at (theta_1=2*pi, theta_2) and (theta_1, theta_2=2*pi)
        Z_i_interp(:,end,1:end-1) = Z_i(:,1,:);     % Add z(theta_1=2*pi, theta_2) = z(theta_1=0, theta_2) (for all theta_2 and for all dimensions of state space)
        Z_i_interp(:,1:end-1,end) = Z_i(:,:,1);     % Add z(theta_1, theta_2=2*pi) = z(theta_1, theta_2=0) (for all theta_1 and for all dimensions of state space)
        Z_i_interp(:,end,end) = Z_i(:,1,1);         % Add z(theta_1=2*pi, theta_2=2*pi) = z(theta_1=0, theta_2=0) (for all dimensions of state space)

        % Interpolate the state vectors on the 2D hypertime grid using csape
        % Here, csape uses periodic end conditions, i.e. the first and second derivatives at the left end are matched with those at the right end (for each grid dimension)
        % To interpolate vectored grid data, the vectors must be stored in a 3-dimensional array with Z_i_interp(:,i,j) being the vector at the grid indices (i,j)
        Z_i_csape = csape({theta_1_interp,theta_2_interp},Z_i_interp,{'periodic','periodic'});
        Z_i_eval = fnval(Z_i_csape,{theta_eval_1,theta_eval_2});        % Evaluate the interpolation at the evaluation coordinates {theta_eval_1,theta_eval_2}

        % Hand the interpolated state vectors Z_i_eval to the output array s
        % size(s) = [res_1 x res_2 x dim x n_evals], but size(Z_i_eval) = [dim x res_1 x res_2]
        if dim == 1                                     % Exception for dim = 1: Z_i_eval is a matrix with size(Z_i_eval) = [res_1 x res_2]
            s(:,:,1,i) = Z_i_eval;                      % Z_i_eval is already good to go and does not have to be permuted
        else
            s(:,:,:,i) = permute(Z_i_eval,[2,3,1]);     % Z_i_eval needs to be permuted so that s has the correct structure
        end

    end


    % Get the mu values
    mu = obj.mu(index);                             % It is ensured by S.solget_up_index that the elements of options.index are unambigious 

end