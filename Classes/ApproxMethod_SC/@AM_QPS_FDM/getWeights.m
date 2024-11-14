% This function is a method of the subclass AM_QPS_FDM.
% It returns the weights w_1_(sigma_k) and w_2_(sigma_k) (as a vector in each case) needed to approximate ...
% dz_i_j/dtheta_1 = 1/DeltaTheta_1 * sum_(k=1)^(p_1) ( w_1_(sigma_1_k) * z_(i+sigma_1_k)_j ) and dz_i_j/dtheta_2 = 1/DeltaTheta_2 * sum_(k=1)^(p_2) ( w_2_(sigma_2_k) * z_i_(j+sigma_2_k) ).
% The required input are either the vectors obj.points_1 and obj.points_2, which store all local grid point indices sigma_k, ...
% or the desired discretization scheme and approximation order. It can be mixed as well, e.g. obj.points_1 for dz_i_j/dtheta_1 and a scheme plus approx. order for dz_i_j/dtheta_2.
% The weights of some typical discretization schemes and approximation orders were deposited in this function.
% If the desired approximation does not fit the deposited cases, the weights are calculated by solving a linear equation system.
%
% @obj:  ApproxMethod subclass AM_QPS_FDM object
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


    scheme_names = fieldnames(avail_schemes);   % Stores the available scheme names



    %% Get the weights (and local grid point indices) for approximation of dz/dtheta_1

    scheme_found_1 = false;                     % This variable states if an implemented scheme has been found during later checks


    % If 'opt_approx_method.points_1' has been given by user
    
    if isfield(DYN.opt_approx_method,'points_1')
        sigma_1 = sort(obj.points_1);                       % Make sure that obj.points_1 is sorted in ascending order (required by check below and QPS_FDM_residuum)
        if size(sigma_1,2) == 1; sigma_1 = sigma_1'; end    % Change sigma_1 to be a row vector if it is a column vector (required to build up the matrix V_1 correctly)
        
        % Check if sigma_1 belongs to an already implemented discretization scheme
        for i = 1:length(scheme_names)         
            for approx_order_1 = 1:size(avail_schemes.(scheme_names{i}),1)
                if isequal(sigma_1,avail_schemes.(scheme_names{i}){approx_order_1,1})   % If an available scheme has been found
                    w_1 = avail_schemes.(scheme_names{i}){approx_order_1,2};            % Set the corresponding weights
                    scheme_1 = scheme_names{i};                                         % Set the scheme (overwrite default setting)
                    scheme_found_1 = true;                                              % Set scheme_found_1 to true in order to terminate the for loops
                    break                                                               % Terminate the inner for loop if a scheme has been found
                end
            end
            if scheme_found_1;   break;   end     % Terminate the outer for loop if a scheme has been found
        end
        
        % If an already implemented discretization scheme has not been found: The weights must be calculated (see below)
        if ~scheme_found_1
            scheme_1 = 'user-defined';                                  % The discretization scheme is user-defined
            approx_order_1 = length(sigma_1) - 1;                       % Not sure if this is true for all cases
            if length(sigma_1(sigma_1>0)) == length(sigma_1(sigma_1<0))
                if sigma_1(sigma_1>0) == sort(abs(sigma_1(sigma_1<0)))  % If user supplied a central scheme
                    scheme_1 = 'central (user-defined)';
                    if isempty(sigma_1(sigma_1==0))                     % If user did not include '0' in 'points'
                        approx_order_1 = approx_order_1 + 1;            % approx_order must be increased by one since grid point index 0 is "missing"
                    end
                end
            end
        end


    % If 'opt_approx_method.scheme_1' and/or 'opt_approx_method.approx_order_1' has been given by user
    % (This part of the code is executed when the user did not input specific grid point indices)
    
    else
        scheme_1 = obj.scheme_1;                % This can always be done as a default scheme is set in AM_QPS_FDM
        approx_order_1 = obj.approx_order_1;    % This can always be done as a default value is set in AM_QPS_FDM
        
        % If approx_order is less or equal than 6 (which is the default value): The grid point indices and weights can be picked from avail_schemes
        if approx_order_1 <= 6                                                  
            sigma_1 = avail_schemes.(scheme_1){approx_order_1,1};       % Set the grid point indices
            w_1 = avail_schemes.(scheme_1){approx_order_1,2};           % Set the corresponding weights
            scheme_found_1 = true;                                      % Set scheme_found to true
        
        % If approx_order is greater than 6: Build the vector of the grid point indices and calculate the weights afterwards
        else
            switch scheme_1
                case 'central';     sigma_1 = -approx_order_1/2:approx_order_1/2;
                case 'forward';     sigma_1 = 0:approx_order_1;
                case 'backward';    sigma_1 = -approx_order_1:0;
            end
        end
        
    end


    % Calculate the weights if necessary
    
    if ~scheme_found_1
        
        p_1 = length(sigma_1);                    % Number of grid points used to approximate dz/dtheta_1
        
        % Build up transponse of the Vandermonde matrix V_1 (note: V_1 is already the transponse of the Vandermonde matrix)
        V_1 = zeros(p_1);
        for k = 1:p_1
            V_1(k,:) = sigma_1.^(k-1);
        end
        
        % Build up right hand side of the equation system V_1 * w_1 = d_1
        d_1 = zeros(p_1,1);
        d_1(2) = 1;
        
        % Solve for the weights
        w_1 = V_1\d_1;
        
    end



    %% Get the weights (and local grid point indices) for discretization of dz/dtheta_2

    scheme_found_2 = false;                     % This variable states if an implemented scheme has been found during later checks

    % If 'opt_approx_method.points_2' has been given by user
    
    if isfield(DYN.opt_approx_method,'points_2')
        sigma_2 = sort(obj.points_2);                       % Make sure that obj.points_2 is sorted in ascending order (required by check below and QPS_FDM_residuum)
        if size(sigma_2,2) == 1; sigma_2 = sigma_2'; end    % Change sigma_2 to be a row vector if it is a column vector (required to build up the matrix V_2 correctly)
        
        % Check if sigma_2 belongs to an already implemented discretization scheme
        for i = 1:length(scheme_names)         
            for approx_order_2 = 1:size(avail_schemes.(scheme_names{i}),1)
                if isequal(sigma_2,avail_schemes.(scheme_names{i}){approx_order_2,1})   % If an available scheme has been found
                    w_2 = avail_schemes.(scheme_names{i}){approx_order_2,2};            % Set the corresponding weights
                    scheme_2 = scheme_names{i};                                         % Set the scheme (overwrite default setting)
                    scheme_found_2 = true;                                              % Set scheme_found_2 to true in order to terminate the for loops
                    break                                                               % Terminate the inner for loop if a scheme has been found
                end
            end
            if scheme_found_2;   break;   end     % Terminate the outer for loop if a scheme has been found
        end
        
        % If an already implemented discretization scheme has not been found: The weights must be calculated (see below)
        if ~scheme_found_2
            scheme_2 = 'user-defined';                                  % The discretization scheme is user-defined
            approx_order_2 = length(sigma_2) - 1;                       % Not sure if this is true for all cases
            if length(sigma_2(sigma_2>0)) == length(sigma_2(sigma_2<0))
                if sigma_2(sigma_2>0) == sort(abs(sigma_2(sigma_2<0)))  % If user supplied a central scheme
                    scheme_2 = 'central (user-defined)';
                    if isempty(sigma_2(sigma_2==0))                     % If user did not include '0' in 'points'
                        approx_order_2 = approx_order_2 + 1;            % approx_order must be increased by one since grid point index 0 is "missing"
                    end
                end
            end
        end


    % If 'opt_approx_method.scheme_2' and/or 'opt_approx_method.approx_order_2' has been given by user
    % (This part of the code is executed when the user did not input specific grid point indices)
    
    else
        scheme_2 = obj.scheme_2;                % This can always be done as a default scheme is set in AM_QPS_FDM
        approx_order_2 = obj.approx_order_2;    % This can always be done as a default value is set in AM_QPS_FDM
        
        % If approx_order is less or equal than 6 (which is the default value): The grid point indices and weights can be picked from avail_schemes
        if approx_order_2 <= 6                                                  
            sigma_2 = avail_schemes.(scheme_2){approx_order_2,1};       % Set the grid point indices
            w_2 = avail_schemes.(scheme_2){approx_order_2,2};           % Set the corresponding weights
            scheme_found_2 = true;                                      % Set scheme_found to true
        
        % If approx_order is greater than 6: Build the vector of the grid point indices and calculate the weights afterwards
        else
            switch scheme_2
                case 'central';     sigma_2 = -approx_order_2/2:approx_order_2/2;
                case 'forward';     sigma_2 = 0:approx_order_2;
                case 'backward';    sigma_2 = -approx_order_2:0;
            end
        end
        
    end


    % Calculate the weights if necessary
    
    if ~scheme_found_2
        
        p_2 = length(sigma_2);                    % Number of grid points used to approximate dz/dtheta_2
        
        % Build up transponse of the Vandermonde matrix V_2 (note: V_2 is already the transponse of the Vandermonde matrix)
        V_2 = zeros(p_2);
        for k = 1:p_2
            V_2(k,:) = sigma_2.^(k-1);
        end
        
        % Build up right hand side of the equation system V_2 * w_2 = d_2
        d_2 = zeros(p_2,1);
        d_2(2) = 1;
        
        % Solve for the weights
        w_2 = V_2\d_2;
        
    end



    %% For the Jacobian matrix
    % Calculate w_1_mat_J = d ( sum_(k=1)^(p_1) ( w_1_(sigma_1_k) * z_(i+sigma_1_k)_j ) ) / d z_i_j  (for all i and for all j) and
    % calculate w_2_mat_J = d ( sum_(k=1)^(p_2) ( w_2_(sigma_2_k) * z_i_(j+sigma_2_k) ) ) / d z_i_j  (for all i and for all j)
    % The matrices w_1_mat_J and w_2_mat_J are two (out of three) components of the derivation d res / ds, which is the major part of the Jacobian matrix
    % This is done here because w_1_mat_J and w_2_mat_J only depend on the weights w_1 and w_2 and are therefore constant throughout the continuation
    % Moreover, the indices p_ind_blkdiag_mat of a block diagonal matrix are calculated here. They are used by the sparse() function (see below)
    
    n_int_1 = obj.n_int_1;                      % Number of hyper-time intervals DeltaTheta_1 into which the hyper-time period 2*pi is divided in theta_1-direction
    n_int_2 = obj.n_int_2;                      % Number of hyper-time intervals DeltaTheta_2 into which the hyper-time period 2*pi is divided in theta_2-direction
    dim = DYN.system.dim;                       % Dimension of state space
    sigma_1_g0 = sigma_1(sigma_1 > 0);          % Get all the elements of sigma_1 which are greater than 0
    sigma_1_s0 = sigma_1(sigma_1 < 0);          % Get all the elements of sigma_1 which are less than 0
    n_sigma_1_g0 = length(sigma_1_g0);          % Get the number of elements of sigma_1 that are greater than 0
    n_sigma_1_s0 = length(sigma_1_s0);          % Get the number of elements of sigma_1 that are less than 0
    sigma_2_g0 = sigma_2(sigma_2 > 0);          % Get all the elements of sigma_2 which are greater than 0
    sigma_2_s0 = sigma_2(sigma_2 < 0);          % Get all the elements of sigma_2 which are less than 0
    n_sigma_2_g0 = length(sigma_2_g0);          % Get the number of elements of sigma_2 that are greater than 0
    n_sigma_2_s0 = length(sigma_2_s0);          % Get the number of elements of sigma_2 that are less than 0
    

    %                 w_1_mat_J                 %
    % w_1_mat_J is a sparse matrix consisting of submatrices w_1_mat located on the main diagonal of w_1_mat_J
    % w_1_mat is equal to the matrix p_w_mat_J which is needed for the Jacobian matrix when calculating periodic solutions (see AM_PS_FDM and the corresponding method getWeights)
    % Therefore, a Matrix W_1 is assembled. The columns of W_1 will be placed on particular diagonals of the matrix w_1_mat
    W_1 = repmat([w_1((end-n_sigma_1_g0+1):end)', w_1', w_1(1:n_sigma_1_s0)'], n_int_1*dim, 1);    
    
    % The next three vectors describe the numbers of the diagonals that are filled by the columns of W_1 (for the numbering, see Matlab doc of spdiags)
    upper_diag_numbers_1 = dim*n_int_1 .* ones(1,n_sigma_1_s0) + dim .* sigma_1_s0;     % This vector sets the diagonal numbers in the upper right corner
    lower_diag_numbers_1 = -dim*n_int_1 .* ones(1,n_sigma_1_g0) + dim .* sigma_1_g0;    % This vector sets the diagonal numbers in the bottom left corner
    diag_numbers_1 = [lower_diag_numbers_1, dim.*sigma_1, upper_diag_numbers_1];        % This vector sets all diagonal numbers  
    
    % Build w_1_mat_J
    w_1_mat = spdiags(W_1, diag_numbers_1, n_int_1*dim, n_int_1*dim);                   % Build w_1_mat using spdiags
    obj.p_w_1_mat_J = kron(eye(n_int_2), w_1_mat);                                      % Place w_1_mat n_int_2-times on the "main diagonal" of w_1_mat_J


    %                 w_2_mat_J                 %
    % w_2_mat_J is a diagonal sparse matrix consisting of the weights w_2 located on certain diagonals
    % w_2_mat_J shares similarities with the Jacobian matrix of periodic solutions w_mat_J
    % In contrast to w_mat_J, w_2_mat_J exhibits (n_int_1-1) zeros(dim)-matrices between the w_2_k.*eye(dim)-matrices when horizontally going though the matrix
    W_2 = repmat([w_2((end-n_sigma_2_g0+1):end)', w_2', w_2(1:n_sigma_2_s0)'], n_int_1*n_int_2*dim, 1); % The columns of W_2 will be placed on particular diagonals of w_mat_2_J

    % The next three vectors describe the numbers of the diagonals that are filled by the columns of W_2 (for the numbering, see Matlab doc of spdiags)
    upper_diag_numbers_2 = dim*n_int_1*n_int_2 .* ones(1,n_sigma_2_s0) + dim*n_int_1 .* sigma_2_s0;     % This vector sets the diagonal numbers in the upper right corner
    lower_diag_numbers_2 = -dim*n_int_1*n_int_2 .* ones(1,n_sigma_2_g0) + dim*n_int_1 .* sigma_2_g0;    % This vector sets the diagonal numbers in the bottom left corner
    diag_numbers_2 = [lower_diag_numbers_2, dim*n_int_1.*sigma_2, upper_diag_numbers_2];                % This vector sets all diagonal numbers 

    % Build w_2_mat using spdiags
    obj.p_w_2_mat_J = spdiags(W_2, diag_numbers_2, n_int_1*n_int_2*dim, n_int_1*n_int_2*dim);


    %              p_ind_blkdiag_mat            %
    % One part to create the the Jacobian is a block diagonal matrix. In order to place the elements at the correct position, the function sparse() is used
    % sparse() needs the row and column indices of the elements which need to be set. These indices are created here
    blkdiag_mat = sparse(kron(eye(n_int_1*n_int_2), ones(dim)));                                                % Create the [n_int_1*n_int_2*dim x n_int_1*n_int_2*dim] block diagonal matrix
    [obj.p_ind_blkdiag_mat(:,1),obj.p_ind_blkdiag_mat(:,2)] = ind2sub(size(blkdiag_mat),find(blkdiag_mat));     % Get the indices of all elements ~= 0 (i.e. = 1)


    
    %% Update important properties and display information

    obj.points_1 = sigma_1;
    obj.points_2 = sigma_2;
    obj.weights_1 = w_1;    
    obj.weights_2 = w_2;

    %{
    disp('Using finite differences to approximate dz(theta_1,theta_2)/dtheta_1 and dz(theta_1,theta_2)/dtheta_2.')
    disp(' ')
    disp(table(categorical({scheme_1;scheme_2}), [approx_order_1;approx_order_2], ...
         'VariableNames', {'Discretization scheme','Approximation order'}, 'RowNames', {'theta_1 direction','theta_2 direction'}))
    disp('The local grid point indices sigma_i_k and corresponding weighting factors w_i_(sigma_k) to approximate dz/dtheta_i (i = 1,2) are:')
    disp(' ')
    disp(table(sigma_1', w_1, 'VariableNames', {'sigma_1_k', 'w_1_(sigma_k)'}))
    disp(table(sigma_2', w_2, 'VariableNames', {'sigma_2_k', 'w_2_(sigma_k)'}))
    disp('-------------------------------------------------------')
    disp('Starting to find initial solution ...')
    %}
    
end