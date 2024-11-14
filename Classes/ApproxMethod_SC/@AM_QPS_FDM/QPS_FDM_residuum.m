% This function is a method of Approx Method subclass AM_QPS_FDM.
% It computes the residuum of the equation system using finite-difference method.
%
% @obj:   ApproxMethod subclass AM_QPS_FDM object
% @y:     Curve point (solution vector to be solved for by fsolve)
% @DYN:   DynamicalSystem class object
%
% @res:   residuum vector of evaluated ODE
% @J_res: Jacobian matrix of res

function [res,J_res] = QPS_FDM_residuum(obj,y,DYN) 

    %% Parameters
    Fcn = DYN.rhs;                              % Right-hand side of the system dz/dt = Fcn(t,z,param)
    dim = DYN.system.dim;                       % Dimension of state space
    n_auto = DYN.n_auto;                        % Number of autonomous frequencies
    
    s = y(1:(end-1-n_auto));                    % s is the FDM solution vector consisting of all discretised state space vectors z(theta_1_i,theta_2_j) (i,j = 0,1,...,n_int_1/2-1)
    mu = y(end);                                % Continuation parameter is the last element of y
    
    if n_auto == 0                              % Non-autonomous system
        Omega = DYN.non_auto_freq(mu);          % Get the angluar frequencies (Omega can be a row vector or a column vector!)
        omega_1 = Omega(1);                     % First non-autonomous frequency
        omega_2 = Omega(2);                     % Second non-autonomous frequency
    elseif n_auto == 1                          % Half-autonomous (mixed) system
        omega_1 = DYN.non_auto_freq(mu);        % Non-autonomous angluar frequency
        omega_2 = y(end-1);                     % Autonomous angluar frequency
        s_p = obj.iv(1:end-1);                  % Solution vector s at predictor point (obj.iv(end) is omega_2 at predictor point)
    elseif n_auto == 2                          % Autonomous system
        omega_1 = y(end-2);                     % First autonomous angluar frequency 
        omega_2 = y(end-1);                     % Second autonomous angluar frequency        
        s_p = obj.iv(1:end-2);                  % Solution vector s at predictor point (obj.iv(end-1:end) is [omega_1; omega_2] at predictor point)
    end
    
    n_int_1 = obj.n_int_1;                      % Number of hyper-time intervals DeltaTheta_1 into which the hyper-time period 2*pi is divided in theta_1-direction
    DeltaTheta_1 = 2*pi / n_int_1;              % Hyper-time interval between two consecutive grid points in theta_1-direction ( DeltaT_1 = DeltaTheta_1 / omega_1 )
    n_int_2 = obj.n_int_2;                      % Number of hyper-time intervals DeltaTheta_2 into which the hyper-time period 2*pi is divided in theta_2-direction
    DeltaTheta_2 = 2*pi / n_int_2;              % Hyper-time interval between two consecutive grid points in theta_2-direction ( DeltaT_2 = DeltaTheta_2 / omega_2 )
    
    sigma_1 = obj.points_1;                     % Local grid point indices sigma_1_k which are needed to approximate dz(theta_1_i,theta_2_j)/dtheta_1
    w_1 = obj.weights_1;                        % Weighting factors which are needed to approximate approximate dz(theta_1_i,theta_2_j)/dtheta_1
    sigma_2 = obj.points_2;                     % Local grid point indices sigma_2_k which are needed to approximate dz(theta_1_i,theta_2_j)/dtheta_2
    w_2 = obj.weights_2;                        % Weighting factors which are needed to approximate approximate dz(theta_1_i,theta_2_j)/dtheta_2

    % Update the continuation parameter in the system parameter array
    param = DYN.param;
    param{DYN.act_param} = mu;

    
    %% Reshape the method solution vector s to a matrix S, assemble an expanded method solution vector matrix S_expnd (simplifies building the equation system) and evaluate f(t,z,param)
    % S = [z(0,0), ... , z(0,theta_2(end)); ...... ; z(theta_1(end),0), ... , z(theta_1(end),theta_2(end))], ...
    % -> the j-th column consists of the state space vectors at (theta_1_i,theta_2_j) with i = 0,1,...,n_int_1 and j = const.
    % S is then expanded to S_expnd = [  zeros() , s_1_before,  zeros() ;
    %                                  s_2_before,      S    , s_2_after;
    %                                    zeros() , s_1_after ,  zeros() ;
    % s_1_before stores all z_i_j that must be placed ahead of z(0,theta_2_j) (sigma_1_k < 0). The z_i_j are taken from the bottom end of S, i.e. z_i_j = z(2*pi + sigma_1_k * DeltaTheta_1, j * DeltaTheta_2)
    % s_1_after stores all z_i_j that must be placed following z(2*pi-DeltaTheta_1,theta_2_j) (sigma_1_k > 0). The z_i_j are taken from the top end of S, i.e. z_i_j = z(sigma_1_k * DeltaTheta_1, j * DeltaTheta_2)
    % s_2_before stores all z_i_j that must be placed ahead of z(theta_1_i,0) (sigma_2_k < 0). The z_i_j are taken from the right end of S, i.e. z_i_j = z(i * DeltaTheta_1, 2*pi + sigma_2_k * DeltaTheta_2)
    % s_2_after stores all z_i_j that must be placed following z(theta_1_i,2*pi-DeltaTheta_2) (sigma_2_k > 0). The z_i_j are taken from the left end of S, i.e. z_i_j = z(i * DeltaTheta_1, sigma_2_k * DeltaTheta_2) 
    % That way, S_expnd already holds the periodic condition in theta_1-direction as well as in theta_2-direction due to its structure
    % This simplifies building the equation system since only a small for-loop is required to get the state space vectors which are needed for building the equation system

    % Reshape s to matrix S
    S = reshape(s,n_int_1*dim,n_int_2);                     % Going down the columns of S represents the theta_1-direction, whereas the theta_2-direction is represented by going from left to right within each row

    % Build S_expnd
    s_1_before = S(end+min(sigma_1)*dim+1 : end, :);        % Get the last abs(min(sigma_1))*dim entries of S in theta_1-direction (theta_1-direction: going down the columns of S)
    s_1_after = S(1 : max(sigma_1)*dim, :);                 % Get the first max(sigma_1)*dim entries of S in theta_1-direction     (theta_1-direction: going down the columns of S)
    s_2_before = S(:, end+min(sigma_2)+1 : end);            % Get the last abs(min(sigma_2))*dim entries of S in theta_2-direction (theta_2-direction: going from left to right within the rows of S)
    s_2_after = S(:, 1 : max(sigma_2));                     % Get the first max(sigma_2)*dim entries of S in theta_2-direction     (theta_2-direction: going from left to right within the rows of S)
    S_expnd = [zeros(size(s_1_before,1),size(s_2_before,2)), s_1_before, zeros(size(s_1_before,1),size(s_2_after,2));       % Assemble S_expnd
                             s_2_before                    ,     S     ,                  s_2_after                 ;       % Assemble S_expnd
               zeros(size(s_1_after,1),size(s_2_before,2)) , s_1_after , zeros(size(s_1_after,1),size(s_2_after,2))];       % Assemble S_expnd

    % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) for every grid point -> theta_1_i = 0, ..., 2*pi-DeltaTheta_1; theta_2_j = 0, ..., 2*pi-DeltaTheta_2 
    % This saves time as no loop is needed for the evaluation of Fcn
    theta_1 = 0 : DeltaTheta_1 : (2*pi-DeltaTheta_1);                   % Create the discretised hyper-time vector theta_1  ->  length(theta_1) = n_int_1
    Theta_1 = repmat(theta_1,1,n_int_2);                                % Theta_1 = [theta_1, ..., theta_1] <n_int_2-times> ->  length(Theta_1) = n_int_1 * n_int_2
    theta_2 = 0 : DeltaTheta_2 : (2*pi-DeltaTheta_2);                   % Create the discretised hyper-time vector theta_2  ->  length(theta_2) = n_int_2
    Theta_2 = reshape(repmat(theta_2,n_int_1,1),1,n_int_1*n_int_2);     % Theta_2 = [theta_2(1) <n_int_1-times>, ... , theta_2(end) <n_int_1-times>]  ->  length(Theta_2) = n_int_1 * n_int_2
    t = [1/omega_1.*Theta_1; 1/omega_2.*Theta_2];                       % Create the discretised time vector for the evaluation of Fcn
    Z = reshape(s,dim,n_int_1*n_int_2);                                 % s is reshaped to Z = [z(0,0), ... , z(theta_1(end),0), z(0,DeltaTheta_2), ... , z(theta_1(end),DeltaTheta_2), ...... , z(theta_1(end),theta_2(end))]
    Fcn_eval = Fcn(t,Z,param);                                          % Fcn_eval is a (dim x n_int_1*n_int_2) array where the (i+1)-th column equals Fcn(t(:,i),Z(:,i),param)


    %% Build up the equation system and the residuum
    % The residuum is res = g in the non-autonomous case where g is the residuum of the approximated ODE:
    % g_i_j = omega_1 .* dz(theta_1_i,theta_2_j)/dtheta_1 + omega_2 .* dz(theta_1_i,theta_2_j)/dtheta_2 - f(t_i_j,z_i_j,param) 
    %       = omega_1/DeltaTheta_1 .* sum_(k=1)^(p_1) ( w_1_(sigma_1_k) .* z_(i+sigma_1_k)_j ) + omega_2/DeltaTheta_2 .* sum_(k=1)^(p_2) ( w_2_(sigma_2_k) .* z_i_(j+sigma_2_k) ) - f(t_i_j,z_i_j,param)
    % for all i = 0:(n_int_1-1), j = 0:(n_int_2-1). This results in n_int_1*n_int_2 equations of dimension dim, leading to length(res) = n_int_1*n_int_2*dim
    % The residuum of the approximated ODE is implemented as: g = omega_1/DeltaTheta_1 .* Z_i_sigma_1 * w_1 + omega_2/DeltaTheta_2 .* Z_j_sigma_2 * w_2 - reshape(Fcn_eval,n_int_1*n_int_2*dim,1)
    % If the system is (half-)autonomous, one or two phase conditions pc2 (and pc1) must be added to g to complete the residuum: res = [g; pc1; pc2]

    Z_i_sigma_1 = zeros(n_int_1*n_int_2*dim,length(sigma_1));           % Preallocate the matrix Z_i_sigma_1 which will store the state space vectors needed to approximate all dz(theta_1_i,theta_2_j)/dtheta_1
    Z_j_sigma_2 = zeros(n_int_1*n_int_2*dim,length(sigma_2));           % Preallocate the matrix Z_j_sigma_2 which will store the state space vectors needed to approximate all dz(theta_1_i,theta_2_j)/dtheta_2 

    idx_shift_1 = size(s_1_before,1)/dim;                               % This is used for an index shift in theta_1-direction
    idx_shift_2 = size(s_2_before,2);                                   % This is used for an index shift in theta_2-direction

    % The derivations dz/dtheta are approximated by dz(theta_1_i,theta_2_j)/dtheta_1 = 1/DeltaTheta_1 .* sum_(k=1)^(p_1) ( w_1_(sigma_1_k) .* z_(i+sigma_1_k)_j ) = 1/DeltaTheta_1 .* Z_i_sigma_1 * w_1 ...
    %                                           and dz(theta_1_i,theta_2_j)/dtheta_2 = 1/DeltaTheta_2 .* sum_(k=1)^(p_2) ( w_2_(sigma_2_k) .* z_i_(j+sigma_2_k) ) = 1/DeltaTheta_2 .* Z_j_sigma_2 * w_2
    % Z_i_sigma_1 and Z_j_sigma_2 are matrices storing all the z_(i+sigma_1_k)_j and z_i_(j+sigma_2_k)
    for k1 = 1:length(sigma_1)
        Z_i_sigma_1(:,k1) = reshape(S_expnd((idx_shift_1+sigma_1(k1))*dim+1 : (idx_shift_1+sigma_1(k1)+n_int_1)*dim, idx_shift_2+1 : idx_shift_2+n_int_2), [n_int_1*n_int_2*dim,1]);   % The columns can be filled by taking submatrices from S_expnd
    end
    for k2 = 1:length(sigma_2)
        Z_j_sigma_2(:,k2) = reshape(S_expnd(idx_shift_1*dim+1 : (idx_shift_1+n_int_1)*dim, (idx_shift_2+sigma_2(k2)+1) : (idx_shift_2+sigma_2(k2)+n_int_2)), [n_int_1*n_int_2*dim,1]);
    end

    % Calculate the residuum of the ODE
    g = omega_1/DeltaTheta_1 .* Z_i_sigma_1 * w_1 + omega_2/DeltaTheta_2 .* Z_j_sigma_2 * w_2 - reshape(Fcn_eval,n_int_1*n_int_2*dim,1);

    % Assemble the residuum res
    if n_auto == 0                                                      % Non-autonomous system
        res = g;                                                        % Build the residuum
    elseif n_auto == 1                                                  % Half-autonomous system 
        pc2 = (Z_j_sigma_2 * w_2)' * (s - s_p);                         % Integral phase condition for autonomous frequency (it is called pc2 here because omega_2 is the autonomous frequency)
        res = [g; pc2];                                                 % Build the residuum
    elseif n_auto == 2                                                  % Full autonomous system 
        pc1 = (Z_i_sigma_1 * w_1)' * (s - s_p);                         % Integral phase condition for first autonomous frequency
        pc2 = (Z_j_sigma_2 * w_2)' * (s - s_p);                         % Integral phase condition for second autonomous frequency
        res = [g; pc1; pc2];                                            % Build the residuum
    end



    %% Calculate the Jacobian matrix J_res
    % Depending on the degree of autonomy, three different Jacobian matrices of res can be built
    % For example, the Jacobian matrix of res in the non-autonomous case consists of J_res = [dg/ds, dg/dmu]
    % In the (half)-autonomous case, the Jacobian matrix features additional derivatives (for the detailed structure, see section "Build the Jacobian matrix J_res")
    % In the following, the derivatives dg/ds and dg/dmu are calculated as they are needed in any case

    h = sqrt(eps);                              % Set the step size used to calculate particular derivatives using forward finite difference
    % h = eps^(1/3);                            % OPTIONAL: Set the step size used to calculate particular derivatives using central finite difference

    % Calculate dg/ds
    % dg_ds = omega_1/DeltaTheta_1 .* d(Z_i_sigma_1 * w_1)/ds + omega_2/DeltaTheta_2 .* d(Z_j_sigma_2 * w_2)/ds - d(reshape(Fcn_eval,n_int_1*n_int_2*dim,1))/ds consists of three main parts: 
    % Two of the three main parts d(Z_i_sigma_1 * w_1)/ds and d(Z_j_sigma_2 * w_2)/ds have already been calculated in getWeights and has been stored in obj.p_w_1_mat_J and obj.p_w_2_mat_J
    % The third main part dFcn_ds_mat = d(reshape(Fcn_eval,n_int_1*n_int_2*dim,1))/ds corresponds to the derivation dFcn(t_i_j)/ds for all i = 0:(n_int_1-1) and j = 0:(n_int_2-1)
    % dFcn_ds_mat is a (n_int_1*n_int_2*dim) x (n_int_1*n_int_2*dim) block diagonal matrix consisting of the derivatives dFcn(t_i_j,z_i_j,param)/dz_i_j placed along its "main diagonal"
    % dFcn_ds_mat is created via the function sparse(). The row and column numbers of the elements to be filled are stored in obj.p_ind_blkdiag_mat, which is created in getWeights
    % The elements of dFcn_ds_mat that are ~= 0 are calculated using forward (OPTIONAL: central) finite difference and provisionally stored in dFcn_ds 
    % The step width to calculate the k-th column of dFcn(t_i_j,z_i_j,param)/dz_i_j is h_(i_j,k) = h*(1+abs(z_(i_j,k)), where z_(i_j,k) is the k-th component of z_i_j
    % In order to save computing time, the code is vectorized and therefore only one (optional: two) evaluation of Fcn is needed
    t_dim = reshape(repmat(t,dim,1),2,n_int_1*n_int_2*dim);             % This is a (2 x n_int_1*n_int_2*dim) vector where each t(:,i) is repeated dim times
    Z_dim = reshape(repmat(Z,dim,1),dim,n_int_1*n_int_2*dim);           % This is a (dim x n_int_1*n_int_2*dim) matrix where each z_i_j is repeated dim times
    H = h.*(repmat(eye(dim),1,n_int_1*n_int_2) + sparse(repmat(1:1:dim,1,n_int_1*n_int_2),1:1:n_int_1*n_int_2*dim,abs(s),dim,n_int_1*n_int_2*dim));     % H stores the individual step widths h_(i_j,k) = h*(1+abs(z_(i_j,k))
    Z_dim_plus_h = Z_dim + H;                                           % Pertub Z_dim by + H
    % Z_dim_minus_h = Z_dim - H;                                        % Pertub Z_dim by - H (OPTIONAL: needed for central finite difference)
    dFcn_ds = (Fcn(t_dim,Z_dim_plus_h,param) - reshape(repmat(Fcn_eval,dim,1),dim,n_int_1*n_int_2*dim)) ./ repmat(nonzeros(H)',dim,1);      % Forward finite difference
    % dFcn_ds = (Fcn(t_dim,Z_dim_plus_h,param) - Fcn(t_dim,Z_dim_minus_h,param)) ./ (2.*repmat(nonzeros(H)',dim,1));                        % OPTIONAL: central finite difference
    dFcn_ds_mat = sparse(obj.p_ind_blkdiag_mat(:,1), obj.p_ind_blkdiag_mat(:,2), reshape(dFcn_ds,n_int_1*n_int_2*dim*dim,1), n_int_1*n_int_2*dim, n_int_1*n_int_2*dim); 
    dg_ds = omega_1/DeltaTheta_1 .* obj.p_w_1_mat_J +  omega_2/DeltaTheta_2 .* obj.p_w_2_mat_J - dFcn_ds_mat;       % Calculate dg/ds 

    % Calculate dg/dmu using forward (OPTIONAL: central) finite difference
    h_mu = h*(1+abs(mu));                                               % Set h_mu to approximate derivatives with respect to mu
    param_plus_h = param;                                               % Define a new parameter array
    param_plus_h{DYN.act_param} = mu + h_mu;                            % Update the new parameter array by mu + h_mu
    % param_minus_h = param;                                            % Define a new parameter array (OPTIONAL: needed for central finite difference)
    % param_minus_h{DYN.act_param} = mu - h_mu;                         % Update the new parameter array by mu - h_mu (OPTIONAL: needed for central finite difference)
    if n_auto == 0                                                      % Non-Autonomous case: mu can be the frequency -> angular frequency omega must be updated
        Omega_plus_h = DYN.non_auto_freq(mu+h_mu);                      % Update the frequencies by "+ h_mu"
        omega_1_plus_h = Omega_plus_h(1);                               % First by "+ h_mu" updated non-autonomous frequency
        omega_2_plus_h = Omega_plus_h(2);                               % Second by "+ h_mu" updated non-autonomous frequency
        % Omega_minus_h = DYN.non_auto_freq(mu-h_mu);                   % Update the frequencies by "- h_mu" (OPTIONAL: needed for central finite difference)
        % omega_1_minus_h = Omega_minus_h(1);                           % First by "- h_mu" updated non-autonomous frequency
        % omega_2_minus_h = Omega_minus_h(2);                           % Second by "- h_mu" updated non-autonomous frequency
        t_plus_h = [1/omega_1_plus_h.*Theta_1; 1/omega_2_plus_h.*Theta_2];          % Update the time vector by "+ h_mu" for the evaluation of Fcn
        % t_minus_h = [1/omega_1_minus_h.*Theta_1; 1/omega_2_minus_h.*Theta_2];     % Update the time vector by "- h_mu" for the evaluation of Fcn (OPTIONAL: needed for central finite difference)
        Fcn_eval_plus_h = Fcn(t_plus_h,Z,param_plus_h);                 % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) with the "mu+h_mu" updated t and param array   
        % Fcn_eval_minus_h = Fcn(t_minus_h,Z,param_minus_h);            % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) with the "mu-h_mu" updated t and param array (OPTIONAL: needed for central finite difference)
        dg_dmu = ( (omega_1_plus_h/DeltaTheta_1.*Z_i_sigma_1*w_1 + omega_2_plus_h/DeltaTheta_2.*Z_j_sigma_2*w_2 - reshape(Fcn_eval_plus_h,n_int_1*n_int_2*dim,1)) - res) / h_mu;    % Forward finite difference 
        % dg_dmu = ( (omega_1_plus_h-omega_1_minus_h)/DeltaTheta_1.*Z_i_sigma_1*w_1 + (omega_2_plus_h-omega_2_minus_h)/DeltaTheta_2.*Z_j_sigma_2*w_2 ...                            % OPTIONAL: central finite difference
        %              - reshape(Fcn_eval_plus_h,n_int_1*n_int_2*dim,1) + reshape(Fcn_eval_minus_h,n_int_1*n_int_2*dim,1) ) / (2*h_mu);                                             % OPTIONAL: central finite difference
    elseif n_auto == 1                                                  % Mixed case: omega_1 is the non-autonomous frequency and therefore it can be a function of mu. omega_2 is the autonomous frequency that is not a function of mu
        omega_1_plus_h = DYN.non_auto_freq(mu+h_mu);                    % Update the non-autonomous frequency by "+ h_mu"
        % omega_1_minus_h = DYN.non_auto_freq(mu-h_mu);                 % Update the non-autonomous frequency by "- h_mu" (OPTIONAL: needed for central finite difference)
        t_plus_h = [1/omega_1_plus_h.*Theta_1; 1/omega_2.*Theta_2];     % Update the time vector by "+ h_mu" for the evaluation of Fcn
        % t_minus_h = [1/omega_1_minus_h.*Theta_1; 1/omega_2.*Theta_2]; % Update the time vector by "- h_mu" for the evaluation of Fcn (OPTIONAL: needed for central finite difference)
        Fcn_eval_plus_h = Fcn(t_plus_h,Z,param_plus_h);                 % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) with the "mu+h_mu" updated t and param array
        % Fcn_eval_minus_h = Fcn(t_minus_h,Z,param_minus_h);            % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) with the "mu-h_mu" updated t and param array (OPTIONAL: needed for central finite difference)
        dg_dmu = ( (omega_1_plus_h/DeltaTheta_1.*Z_i_sigma_1 * w_1 + omega_2/DeltaTheta_2.*Z_j_sigma_2 * w_2 - reshape(Fcn_eval_plus_h,n_int_1*n_int_2*dim,1)) - res(1:n_int_1*n_int_2*dim)) / h_mu;    % Forward finite difference 
        % dg_dmu = ( (omega_1_plus_h-omega_1_minus_h)/DeltaTheta_1.*Z_i_sigma_1 * w_1 - reshape(Fcn_eval_plus_h,n_int_1*n_int_2*dim,1) + reshape(Fcn_eval_minus_h,n_int_1*n_int_2*dim,1) ) / (2*h_mu);  % OPTIONAL: central finite difference
    elseif n_auto == 2                                                  % Full autonomous case: Both omega_1 and omega_2 are not a function of mu
        Fcn_eval_plus_h = Fcn(t,Z,param_plus_h);                        % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) with the "mu+h" updated t and param array
        % Fcn_eval_minus_h = Fcn(t,Z,param_minus_h);                    % Evaluate the rhs of dz/dtheta_1 .* omega_1 + dz/dtheta_2 .* omega_2 = f(t,z,param) with the "mu-h" updated t and param array (OPTIONAL: needed for central finite difference)
        dg_dmu = ( - reshape(Fcn_eval_plus_h,n_int_1*n_int_2*dim,1) + reshape(Fcn_eval,n_int_1*n_int_2*dim,1) ) / h_mu;                   % Forward finite difference 
        % dg_dmu = ( - reshape(Fcn_eval_plus_h,n_int_1*n_int_2*dim,1) + reshape(Fcn_eval_minus_h,n_int_1*n_int_2*dim,1) ) / (2*h_mu);     % OPTIONAL: central finite difference
    end

    
    % Build the Jacobian matrix J_res
    if n_auto == 0                                                      % Non-autonomous system
         J_res = [dg_ds, dg_dmu];                                       % Build the Jacobian matrix
    elseif n_auto == 1                                                  % Half-autonomous system: omega_2 is the autonomous frequency (here, the phase condition is called pc2)
        dg_domega2 = 1/DeltaTheta_2 .* Z_j_sigma_2 * w_2;               % dg/domega2
        dpc2_ds = (s - s_p)' * obj.p_w_2_mat_J + (Z_j_sigma_2 * w_2)';  % dpc2/ds
        dpc2_domega2 = 0;     dpc2_dmu = 0;                             % Phase condition is independent of the autonomous frequency and of the continuation parameter
        J_res = [dg_ds,   dg_domega2,   dg_dmu;                         % Build the Jacobian matrix
                 dpc2_ds, dpc2_domega2, dpc2_dmu];                      % Build the Jacobian matrix
    elseif n_auto == 2                                                  % Full autonomous system
        dg_domega1 = 1/DeltaTheta_1 .* Z_i_sigma_1 * w_1;               % dg/domega1 
        dg_domega2 = 1/DeltaTheta_2 .* Z_j_sigma_2 * w_2;               % dg/domega2
        dpc1_ds = (s - s_p)' * obj.p_w_1_mat_J + (Z_i_sigma_1 * w_1)';  % dpc1/ds
        dpc2_ds = (s - s_p)' * obj.p_w_2_mat_J + (Z_j_sigma_2 * w_2)';  % dpc2/ds
        dpc1_domega1 = 0;     dpc1_domega2 = 0;     dpc1_dmu = 0;       % Phase condition in theta_1-direction is independent of the frequencies and of the continuation parameter
        dpc2_domega1 = 0;     dpc2_domega2 = 0;     dpc2_dmu = 0;       % Phase condition in theta_2-direction is independent of the frequencies and of the continuation parameter         
        J_res = [dg_ds,   dg_domega1,   dg_domega2,   dg_dmu;           % Build the Jacobian matrix
                 dpc1_ds, dpc1_domega1, dpc1_domega2, dpc1_dmu;         % Build the Jacobian matrix
                 dpc2_ds, dpc2_domega1, dpc2_domega2, dpc2_dmu; ];      % Build the Jacobian matrix
    end

end