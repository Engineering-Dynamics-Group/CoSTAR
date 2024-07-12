% This function is a method of Approx Method subclass AM_PS_FDM.
% It computes the residuum of the equation system using finite-difference method.
%
% @obj:   ApproxMethod subclass AM_PS_FDM object
% @y:     Curve point (solution vector to be solved for by fsolve)
% @DYN:   DynamicalSystem class object
%
% @res:   residuum vector of evaluated ODE
% @J_res: Jacobian matrix of res

function [res,J_res] = PS_FDM_residuum(obj,y,DYN) 

    %% Parameters
    dim = DYN.system.dim;                       % Dimension of state space
    Fcn = DYN.rhs;                              % Right-hand side of the system dz/dt = Fcn(t,z,param)
    n_auto = DYN.n_auto;                        % Number of autonomous frequencies
        
    s = y(1:(end-1-n_auto));                    % s is the FDM solution vector consisting of all discretised state space vectors z(theta_i) (i = 0,1,...,n_int-1)
    mu = y(end);                                % Continuation parameter is the last element of y
    
    if n_auto == 0
        omega = DYN.non_auto_freq(mu);          % Angluar frequency (non-autonomous system)  ( periodic time: T = 2*pi / omega )
    elseif n_auto == 1
        omega = y(end-1);                       % Angluar frequency (autonomous system)  ( periodic time: T = 2*pi / omega )
        s_p = obj.iv(1:end-1);                  % Solution vector s at the predictor point (obj.iv(end) is the autonomous frequency at the predictor point)
    end
    
    n_int = obj.n_int;                          % Number of hyper-time intervals DeltaTheta into which the hyper-time period 2*pi is divided
    DeltaTheta = 2*pi / n_int;                  % Distance in hyper-time coordinate theta between two consecutive grid points for dicretization  ( DeltaT = DeltaTheta / omega )
    
    sigma = obj.points;                         % Local grid point indices sigma_k which are needed to approximate dz/dtheta at theta_i
    w = obj.p_weights;                          % Weighting factors needed to approximate dz/dtheta at theta_i

    % Update the continuation parameter in the system parameter array
    param = DYN.param;
    param{DYN.act_param} = mu;

    
    %% Assemble an expanded method solution vector s_expnd (simplifies building the equation system) and evaluate the rhs of dz/dtheta .* omega = f(t,z,param)
    % s_expnd = [s_before; s; s_after];
    % s_before stores all z_i that must be placed ahead of s (sigma_k < 0). The z_i are taken from the bottom end of s, i.e. z_i = z(2*pi + sigma_k * DeltaTheta)
    % s_after stores all z_i that must be placed following s (sigma_k > 0). The z_i are taken from the top end of s, i.e. z_i = z(sigma_k * DeltaTheta)
    % That way, s_expnd already holds the periodic condition due to its structure
    % This simplifies building the equation system since only a small for-loop is required to get the state space vectors which are needed for building the equation system

    % Build s_expnd
    s_before = s(end+min(sigma)*dim+1 : end);       % Get the last abs(min(sigma))*dim entries of s that are necessary to build up the equation system
    s_after = s(1 : max(sigma)*dim);                % Get the first max(sigma)*dim entries of s that are necessary to build up the equation system
    s_expnd = [s_before; s; s_after];               % Assemble s_expnd

    % Evaluate the rhs of dz/dtheta .* omega = f(t,z,param) for all theta_i in [0, 2*pi-DeltaTheta]. This saves time as Fcn does not have to be called in every loop
    theta = 0 : DeltaTheta : 2*pi - DeltaTheta;     % Build the theta vector for the evaluation  ( t = 1/omega.*theta )
    Z = reshape(s,dim,n_int);                       % s is reshaped to [z(0), z(DeltaTheta), ... , z(2*pi-DeltaTheta)]
    Fcn_eval = Fcn(1/omega.*theta,Z,param);         % Fcn_eval is a (dim x n_int) array where the (i+1)-th colums equals F(t_i,z_i,param)


    %% Build up the equation system and the residuum
    % The residuum is res = g in the non-autonomous case and res = [g; pc] in the autonomous case, where g is the residuum of the approximated ODE:
    % g_i = omega .* dz(theta_i)/dtheta - f(t_i,z_i,param) = omega/DeltaTheta .* sum_(k=1)^p ( w_(sigma_k) .* z_(i+sigma_k) ) - f(t_i,z_i,param)
    % for all i = 0:(n_int-1). This results in n_int equations of dimension dim, leading to length(g) = n_int*dim
    % The residuum of the approximated ODE is implemented as: g = omega/DeltaTheta .* Z_i_sigma * w - reshape(Fcn_eval,n_int*dim,1)
    % pc is the phase condition, which is required in the autonomous case
    
    Z_i_sigma = zeros(n_int*dim,length(sigma));                         % Preallocate the matrix Z_i_sigma which will store the state space vectors needed to approximate all dz(theta_i)/dtheta

    idx_shift = length(s_before)/dim;                                   % This is used for an index shift -> idx_shift = number of elements of sigma < 0

    % dz/dtheta is approximated by dz(theta_i)/dt = 1/DeltaTheta .* sum_(k=1)^p ( w_k .* z_(i+sigma_k) ) = 1/DeltaTheta .* Z_i_sigma * w
    % Z_i_sigma is a (n_int*dim x length(sigma)) matrix storing all the z_(i+sigma_k)
    for i = 1:length(sigma)
        Z_i_sigma(:,i) = s_expnd((idx_shift+sigma(i))*dim+1 : (idx_shift+sigma(i)+n_int)*dim);      % The columns can directly be taken from s_expnd 
    end 

    % Calculate the residuum of the ODE
    g = omega/DeltaTheta .* Z_i_sigma * w - reshape(Fcn_eval,n_int*dim,1);        

    % Assemble the residuum res
    if n_auto == 0                                                      % Non-autonomous system
        res = g;                                                        % Build the residuum
    elseif n_auto == 1                                                  % Autonomous system
        pc = (Z_i_sigma * w)' * (s - s_p);                              % Integral phase condition: sum_(i=0)^(n_int-1) ( f(t_i,z_i,param)' * (z(theta_i) - z_p(theta_i)) ), but ...
        % Poincare phase condition (not used anymore)                   % ... f(t_i,z_i,param) = dz(theta_i)/dtheta .* omega is approximated by FD -> benefit: dpc/dmu = 0
        %{
        z_0 = s(1:dim);                                                 % State space vector at theta = 0
        zp_0 = s_p(1:dim);                                              % State space vector at the predictor point for theta = 0
        Fcn_eval_pc = Fcn(0,zp_0,param);                                % Evaluate Fcn for the phase condition
        pc = Fcn_eval_pc' * (z_0 - zp_0);                               % Poincare phase condition
        %}
        res = [g; pc];                                                  % Build the residuum
    end


    %% Calculate the Jacobian matrix J_res
    % The Jacobian matrix of res consists of J_res = [dg/ds, dg/dmu] in the non-autonomous case ...
    % and consists of J_res = [dg/ds, dg/domega, dg/dmu; dpc/ds, dpc/domega, dpc/dmu] in the autonomous case, ...
    % where pc is the phase condition and where omega is the autonomous frequency (for g, s and mu, see the other explanations in this script).
    % In the following, the derivatives dg/ds and dg/dmu are calculated as they are needed in any case
    
    h = sqrt(eps);                              % Set the step size used to calculate particular derivatives using forward finite difference
    % h = eps^(1/3);                            % OPTIONAL: Set the step size used to calculate particular derivatives using central finite difference

    % Calculate dg/ds
    % dg_ds = omega/DeltaTheta .* d(Z_i_sigma * w)/ds - d(reshape(Fcn_eval,n_int*dim,1))/ds consists of two main parts: 
    % The first main part d(Z_i_sigma * w)/ds has already been calculated in getWeights and has been stored in obj.p_w_mat_J
    % The second main part dFcn_ds_mat = d(reshape(Fcn_eval,n_int*dim,1))/ds corresponds to the derivation dFcn(t_i)/ds for all i = 0:(n_int-1) 
    % dFcn_ds_mat is a (n_int*dim) x (n_int*dim) block diagonal matrix consisting of the derivatives dFcn(t_i,z_i,param)/dz_i placed along its "main diagonal"
    % dFcn_ds_mat is created via the function sparse(). The row and column numbers of the elements to be filled are stored in obj.p_ind_blkdiag_mat, which is created in getWeights
    % The elements of dFcn_ds_mat that are ~= 0 are calculated using forward (OPTIONAL: central) finite difference and provisionally stored in dFcn_ds
    % The step width to calculate the k-th column of dFcn(t_i,z_i,param)/dz_i is h_(i,k) = h*(1+abs(z_(i,k)), where z_(i,k) is the k-th component of z_i
    % In order to save computing time, the code is vectorized and therefore only one (OPTIONAL: two) evaluation of Fcn is needed
    theta_dim = reshape(repmat(theta,dim,1),1,n_int*dim);               % This is a (1 x n_int*dim) vector where each theta_i is repeated dim times
    Z_dim = reshape(repmat(Z,dim,1),dim,n_int*dim);                     % This is a (dim x n_int*dim) matrix where each z_i is repeated dim times
    H = h.*(repmat(eye(dim),1,n_int) + sparse(repmat(1:1:dim,1,n_int),1:1:n_int*dim,abs(s),dim,n_int*dim));     % H stores the individual step widths h_(i,j) = h*(1+abs(Z_(i,j))
    Z_dim_plus_h = Z_dim + H;                                           % Pertub Z_dim by + H
    % Z_dim_minus_h = Z_dim - H;                                        % Pertub Z_dim by - H (OPTIONAL: needed for central finite difference)
    dFcn_ds = (Fcn(1/omega.*theta_dim,Z_dim_plus_h,param) - reshape(repmat(Fcn_eval,dim,1),dim,n_int*dim)) ./ repmat(nonzeros(H)',dim,1);           % Forward finite difference
    % dFcn_ds = (Fcn(1/omega.*theta_dim,Z_dim_plus_h,param) - Fcn(1/omega.*theta_dim,Z_dim_minus_h,param)) ./ (2.*repmat(nonzeros(H)',dim,1));      % OPTIONAL: central finite difference
    dFcn_ds_mat = sparse(obj.p_ind_blkdiag_mat(:,1), obj.p_ind_blkdiag_mat(:,2), reshape(dFcn_ds,n_int*dim*dim,1), n_int*dim, n_int*dim); 
    dg_ds = omega/DeltaTheta .* obj.p_w_mat_J - dFcn_ds_mat;            % Calculate dg/ds 

    % Calculate dg/dmu using forward (OPTIONAL: central) finite difference (attention: the phase condition must be excluded from res)
    h_mu = h*(1+abs(mu));                                               % Set h_mu to approximate derivatives with respect to mu
    param_plus_h = param;                                               % Define a new parameter array
    param_plus_h{DYN.act_param} = mu + h_mu;                            % Update the new parameter array by mu + h_mu
    % param_minus_h = param;                                            % Define a new parameter array (OPTIONAL: needed for central finite difference)
    % param_minus_h{DYN.act_param} = mu - h_mu;                         % Update the new parameter array by mu - h_mu (OPTIONAL: needed for central finite difference)
    if n_auto == 0                                                      % Non-Autonomous case: mu can be the frequency -> angular frequency omega must be updated
        omega_plus_h = DYN.non_auto_freq(mu+h_mu);                      % Update the frequency
        % omega_minus_h = DYN.non_auto_freq(mu-h_mu);                   % Update the frequency (OPTIONAL: needed for central finite difference)
        Fcn_eval_plus_h = Fcn(1/omega_plus_h.*theta,Z,param_plus_h);        % Evaluate the rhs of omega.*dz/dtheta = f(t,z,param) with the "mu+h_mu" updated t and param array   
        % Fcn_eval_minus_h = Fcn(1/omega_minus_h.*theta,Z,param_minus_h);   % Evaluate the rhs of omega.*dz/dtheta = f(t,z,param) with the "mu-h_mu" updated t and param array (OPTIONAL: needed for central finite difference)
        dg_dmu = ( (omega_plus_h/DeltaTheta.*Z_i_sigma*w - reshape(Fcn_eval_plus_h,n_int*dim,1)) - g) / h_mu;          % Forward finite difference
        % dg_dmu = ( (omega_plus_h/DeltaTheta.*Z_i_sigma*w - reshape(Fcn_eval_plus_h,n_int*dim,1)) - (omega_minus_h/DeltaTheta.*Z_i_sigma*w - reshape(Fcn_eval_minus_h,n_int*dim,1))) / (2*h_mu);      % OPTIONAL: central finite difference
    elseif n_auto == 1                                                  % Autonomous case: mu is not the frequency and omega does not have to be updated                                               
        Fcn_eval_plus_h = Fcn(1/omega.*theta,Z,param_plus_h);           % Evaluate the rhs of omega.*dz/dtheta = f(t,z,param) with the "mu+h_mu" updated param array (t does not need to be updated since system is autonomous)
        % Fcn_eval_minus_h = Fcn(1/omega.*theta,Z,param_minus_h);       % Evaluate the rhs of omega.*dz/dtheta = f(t,z,param) with the "mu-h_mu" updated t and param array (OPTIONAL: needed for central finite difference)
        dg_dmu = ( - reshape(Fcn_eval_plus_h,n_int*dim,1) + reshape(Fcn_eval,n_int*dim,1) ) / h_mu;                 % Forward finite difference. OPTIONAL: Use h = h*(1+abs(mu))  
        % dg_dmu = ( - reshape(Fcn_eval_plus_h,n_int*dim,1) + reshape(Fcn_eval_minus_h,n_int*dim,1) ) / (2*h_mu);   % OPTIONAL: central finite difference
    end

    % Build the Jacobian matrix J_res
    if n_auto == 0                                                      % Non-autonomous system
        J_res = [dg_ds, dg_dmu];                                        % Build the Jacobian matrix
    elseif n_auto == 1                                                  % Autonomous system: Additional derivations must be calculated   
        dg_domega = 1 / DeltaTheta .* Z_i_sigma * w;                    % dg/domega 
        dpc_ds = (s - s_p)' * obj.p_w_mat_J + (Z_i_sigma * w)';         % dpc/ds (pc: integral phase condition, see above)
        % dpc_ds = [Fcn_eval_pc', zeros(1,(n_int-1)*dim)];              % dpc/ds of Poincare phase condition (not used anymore, see above)
        dpc_domega = 0;                                                 % dpc/domega = 0, because the phase condition is independent of the frequency
        dpc_dmu = 0;                                                    % dpc/dmu = 0, because the phase condition is independent of the continuation parameter
        % dpc_dmu = (Fcn(0,zp_0,param_plus_h)' - Fcn_eval_pc') * (z_0 - zp_0) / h;                      % dpc/dmu of Poincare phase condition calculated using forward finite difference (not used anymore, see above)
        % dpc_dmu = (Fcn(0,zp_0,param_plus_h)' - Fcn(0,zp_0,param_minus_h)') * (z_0 - zp_0) / (2*h);    % OPTIONAL: dpc/dmu of Poincare phase condition calculated using central finite difference (not used anymore, see above)
        J_res = [dg_ds,  dg_domega,  dg_dmu;                            % Build the Jacobian matrix
                 dpc_ds, dpc_domega, dpc_dmu];                          % Build the Jacobian matrix
    end
    

end