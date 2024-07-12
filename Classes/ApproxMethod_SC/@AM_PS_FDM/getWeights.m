% This function is a method of the subclass AM_PS_FDM.
% It returns the weights w_(sigma_k) (as a vector) needed to approximate dz(theta_i)/dtheta = 1/DeltaTheta * sum_(k=1)^p ( w_(sigma_k) * z_(i+sigma_k) ).
% The required input is either the vector obj.points, which stores all indices sigma_k, or the desired discretization scheme and approximation order.
% The weights of some typical discretization schemes and approximation orders were deposited in this function.
% If the desired approximation does not fit the deposited cases, the weights are calculated by solving a linear equation system.
%
% @obj:  ApproxMethod subclass AM_PS_FDM object
% @DYN:  DynamicalSystem class object

function obj = getWeights(obj,DYN)
    
    %% Implemented discretization schemes and approximation orders
    % The grid point indices and weights of common central, forward and backward discretization schemes of orders 1 to 6 are already implemented
    % The first column each stores the grid point indices sigma_k, the second column each stores the corresponding weights w_(sigma_k)
    % Each row number corresponds to the approximation order of that discretization scheme
    % 'central': Uneven approximation orders are not possible, but the code requires the rows to be present, which is why some values are NaN 

    avail_schemes.central = {           NaN,                               NaN                  ;       % not possible
                             [       -1, 0, 1      ],   [             -1/2; 0; 1/2             ];       % approx_order: 2
                                        NaN,                               NaN                  ;       % not possible
                             [   -2, -1, 0, 1, 2   ],   [       1/12; -2/3; 0; 2/3; -1/12      ];       % approx_order: 4
                                        NaN,                               NaN                  ;       % not possible
                             [-3,-2, -1, 0, 1, 2, 3],   [-1/60; 3/20; -3/4; 0; 3/4; -3/20; 1/60]};      % approx_order: 6
    
    avail_schemes.forward = {[0, 1               ],   [     -1;  1                                 ];   % approx_order: 1
                             [0, 1, 2            ],   [   -3/2;  2;  -1/2                          ];   % approx_order: 2
                             [0, 1, 2, 3         ],   [  -11/6;  3;  -3/2;   1/3                   ];   % approx_order: 3
                             [0, 1, 2, 3, 4      ],   [ -25/12;  4;    -3;   4/3;  -1/4            ];   % approx_order: 4
                             [0, 1, 2, 3, 4, 5   ],   [-137/60;  5;    -5;  10/3;  -5/4;  1/5      ];   % approx_order: 5
                             [0, 1, 2, 3, 4, 5, 6],   [ -49/20;  6; -15/2;  20/3; -15/4;  6/5; -1/6]};  % approx_order: 6

    avail_schemes.backward = {[                    -1, 0],   [                                 -1;       1];    % approx_order: 1
                              [                -2, -1, 0],   [                            1/2; -2;     3/2];    % approx_order: 2
                              [            -3, -2, -1, 0],   [                    -1/3;   3/2; -3;    11/6];    % approx_order: 3
                              [        -4, -3, -2, -1, 0],   [              1/4;  -4/3;     3; -4;   25/12];    % approx_order: 4
                              [    -5, -4, -3, -2, -1, 0],   [      -1/5;   5/4; -10/3;     5; -5;  137/60];    % approx_order: 5
                              [-6, -5, -4, -3, -2, -1, 0],   [ 1/6; -6/5;  15/4; -20/3;  15/2; -6;   49/20]};   % approx_order: 6 

    scheme_found = false;                       % This variable states if an implemented scheme has been found during later checks
    

    %% If 'opt_approx_method.points' has been given by user

    if isfield(DYN.opt_approx_method,'points')
        sigma = sort(obj.points);                           % Make sure that obj.points is sorted in ascending order (required by check below and PS_FDM_residuum)
        if size(sigma,2) == 1; sigma = sigma'; end          % Change sigma to be a row vector if it is a column vector (required to build up the matrix V correctly)

        % Check if sigma belongs to an already implemented discretization scheme
        scheme_names = fieldnames(avail_schemes);
        for i = 1:length(scheme_names)         
            for approx_order = 1:size(avail_schemes.(scheme_names{i}),1)
                if isequal(sigma,avail_schemes.(scheme_names{i}){approx_order,1})       % If an available scheme has been found
                    w = avail_schemes.(scheme_names{i}){approx_order,2};                % Set the corresponding weights
                    scheme = scheme_names{i};                                           % Set the scheme (overwrite default setting)
                    scheme_found = true;                                                % Set scheme_found to true in order to terminate the for loops
                    break                                                               % Terminate the inner for loop if a scheme has been found
                end
            end
            if scheme_found;   break;   end     % Terminate the outer for loop if a scheme has been found
        end

        % If an already implemented discretization scheme has not been found: The weights must be calculated (see below)
        if ~scheme_found
            scheme = 'user-defined';                                    % The discretization scheme is user-defined
            approx_order = length(sigma) - 1;                           % Not sure if this is true for all cases
            if length(sigma(sigma>0)) == length(sigma(sigma<0)) 
                if sigma(sigma>0) == sort(abs(sigma(sigma<0)))          % If user supplied a central scheme
                    scheme = 'central (user-defined)';
                    if isempty(sigma(sigma==0))                         % If user did not include '0' in 'points'
                        approx_order = approx_order + 1;                % approx_order must be increased by one since grid point index 0 is "missing"
                    end
                end
            end
        end


    %% If 'opt_approx_method.scheme' and/or 'opt_approx_method.approx_order' has been given by user
    % This part of the code is executed when the user did not input specific grid point indices

    else
        scheme = obj.scheme;                    % This can always be done as a default scheme is set in AM_PS_FDM
        approx_order = obj.approx_order;        % This can always be done as a default value is set in AM_PS_FDM
        
        % If approx_order is less or equal than 6 (which is the default value): The grid point indices and weights can be picked from avail_schemes
        if approx_order <= 6                                                  
            sigma = avail_schemes.(scheme){approx_order,1};             % Set the grid point indices
            w = avail_schemes.(scheme){approx_order,2};                 % Set the corresponding weights
            scheme_found = true;                                        % Set scheme_found to true
        
        % If approx_order is greater than 6: Build the vector of the grid point indices and calculate the weights afterwards
        else
            switch scheme
                case 'central';     sigma = -approx_order/2:approx_order/2;
                case 'forward';     sigma = 0:approx_order;
                case 'backward';    sigma = -approx_order:0;
            end
        end
    
    end


    %% Calculate the weights if necessary

    if ~scheme_found

        p = length(sigma);                     % Number of grid points used to approximate dz/dt

        % Build up transponse of the Vandermonde matrix V (note: V is already the transponse of the Vandermonde matrix)
        V = zeros(p);
        for k = 1:p
            V(k,:) = sigma.^(k-1);
        end

        % Build up right hand side of the equation system V * w = d
        d = zeros(p,1);
        d(2) = 1;

        % Solve for the weights
        w = V\d;

    end


    %% For Jacobian matrix
    % Calculate w_mat_J = d ( sum_(k=1)^p ( w_(sigma_k) * z_(i+sigma_k) ) ) / d z_i  (for all i)
    % The matrix w_mat_J is one component of the derivation d res / ds, which is the major part of the Jacobian matrix
    % This is done here because w_mat_J only depends on the weights w and is therefore constant throughout the continuation
    % Moreover, the indices p_ind_blkdiag_mat of a block diagonal matrix are calculated here. They are used by the sparse() function (see below)
    
    n_int = obj.n_int;                          % Number of intervals in [0, 2*pi] used for the finte-difference discretization
    dim = DYN.system.dim;                       % Dimension of state space
    sigma_g0 = sigma(sigma > 0);
    sigma_s0 = sigma(sigma < 0);
    n_sigma_g0 = length(sigma_g0);              % Get the number of elements of sigma that are greater than 0
    n_sigma_s0 = length(sigma_s0);              % Get the number of elements of sigma that are smaller than 0
    

    %                  w_mat_J                  %
    % w_mat_J is a diagonal sparse matrix with elements near (or on) the main diagonal, in the bottom left corner and in the upper right corner
    % Therefore, a Matrix W is assembled. The columns of W will be placed on particular diagonals of the matrix w_mat_J
    W = repmat([w((end-n_sigma_g0+1):end)', w', w(1:n_sigma_s0)'], n_int*dim, 1);    
    
    % The next three vectors describe the number of the diagonals that are filled by the columns of W (for the numbering, see Matlab doc of spdiags)
    upper_diag_numbers = dim*n_int * ones(1,n_sigma_s0) + dim .* sigma_s0;      % This vector sets the diagonal numbers in the upper right corner
    lower_diag_numbers = -dim*n_int * ones(1,n_sigma_g0) + dim .* sigma_g0;     % This vector sets the diagonal numbers in the bottom left corner
    diag_numbers = [lower_diag_numbers, dim.*sigma, upper_diag_numbers];        % This vector sets all diagonal numbers  
    
    % Finally, obj.p_w_mat_J = w_mat_J is built up using spdiags
    obj.p_w_mat_J = spdiags(W, diag_numbers, n_int*dim, n_int*dim);


    %              p_ind_blkdiag_mat            %
    % One part to create the the Jacobian is a block diagonal matrix. In order to place the elements at the correct position, the function sparse() is used
    % sparse() needs the row and column indices of the elements which need to be set. These indices are created here
    blkdiag_mat = sparse(kron(eye(n_int), ones(dim)));                                                          % Create the [n_int*dim x n_int*dim] block diagonal matrix
    [obj.p_ind_blkdiag_mat(:,1),obj.p_ind_blkdiag_mat(:,2)] = ind2sub(size(blkdiag_mat),find(blkdiag_mat));     % Get the indices of all elements ~= 0 (i.e. = 1)



    %% Update important properties and display information

    obj.points = sigma;
    obj.p_weights = w;

    disp('Using finite differences to approximate dz(theta)/dtheta')
    disp(' ')
    disp(table(categorical({scheme}), approx_order, 'VariableNames', {'Discretization scheme','Approximation order'}))
    % disp(append('Discretization scheme: ', scheme))
    % disp(append('Approximation order: ', num2str(approx_order)))
    disp('The local grid point indices sigma_k and corresponding weighting factors w_(sigma_k) to approximate dz/dtheta are:')
    disp(' ')
    disp(table(sigma', w, 'VariableNames', {'sigma_k', 'w_(sigma_k)'}))
    disp('-------------------------------------------------------')
    disp('Starting to find initial solution ...')
    

end