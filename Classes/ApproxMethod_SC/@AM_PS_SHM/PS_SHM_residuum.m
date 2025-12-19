% Function for generating the residuum by means of periodic multiple shooting method
%
% @obj:   ApproxMethod subclass AM_PS_SHM object
% @y:     solution vector for the continuation (contains continuation parameter)
% @DYN:   DynamicalSystem class object
%
% @res:   residual vector for the newton-type corrector method in Continuation class
% @J_res: Jacobian of res


function [res,J_res] = PS_SHM_residuum(obj,y,DYN)

    % Parameter (defining these makes the code faster)
    dim = DYN.dim;                                      % Dimension of the system
    Fcn = DYN.rhs;                                      % RHS of the system
    n_auto = DYN.n_auto;                                % Number of autonomous frequencies
    n_shoot = obj.n_shoot;                              % Number of shooting points
    n_time = ceil(100/obj.n_shoot);                     % Number of time evaluation points in each shooting interval for the integral phase condition
    s = y(1:(end-1-n_auto));                            % Method solution vector (hyper-vector of shooting points [z_0; z_1; ...; z_n])
    mu = y(end);                                        % Continuation parameter

    if n_auto == 0
        omega = DYN.non_auto_freq(mu);                  % Angular frequency (non-autonomous system)
    elseif n_auto == 1
        omega = y(end-1);                               % Angular frequency (autonomous system)
    end

    param = DYN.param;                                  % Parameter array
    param{DYN.act_param} = mu;                          % Update the active parameter

    calc_stability = strcmpi(DYN.stability,'on');       % Logical variable defining if stability is computed
    if calc_stability;  h = eps^(1/3);                  % When stability is computed: Jacobian is approximated using central finite difference
    else;               h = sqrt(eps);                  % When stability is not computed: Jacobian is approximated using forward finite difference
    end


    % Preparation for integration
    T = 2*pi./omega;                                    % Periodic time
    dT = T/n_shoot;                                     % Time span for each shooting operation (integration)
    T_int = [0:dT:(n_shoot-1)*dT; dT:dT:n_shoot*dT].';  % Define time intervals for the shooting operation (integration)
    z0_mat = reshape(s,dim,n_shoot);                    % Reshape s to a matrix of size [dim x n_shoot]
    s_p_mat = reshape(obj.iv(1:end-n_auto),dim,n_shoot);% Matrix of shooting points at the predictor point (obj.iv(end) is the autonomous frequency if present)
    Z_traj = zeros(dim,n_time,n_shoot);                 % Stores the trajectory for evaluating the integral phase condition
    Z_end = zeros(dim,n_shoot);                         % Stores the end points of integration (could be stored in Z_traj as well, but using Z_end is more convenient)
    Z_p = zeros(dim,n_time,n_shoot);                    % Stores the trajectory of the predictor point for evaluating the integral phase condition
    odeOpts_1 = obj.odeOpts;                            % Use these options for the integrations with perturbed z_i (because that is where the FCN_wrapper is used)
    odeOpts_2 = obj.odeOpts;                            % Use these options for the integrations with perturbed mu (FCN_wrapper is used as well)
    if strcmpi(obj.solver,'ode15s') || strcmpi(obj.solver,'ode23s') || strcmpi(obj.solver,'ode23t') || strcmpi(obj.solver,'ode23tb')
        odeOpts_1.JPattern = kron(speye(2*dim+2),spones(ones(dim)));    % Specify the Jacobian pattern for implicit solvers (used for time step, not corrector step)
        odeOpts_2.JPattern = kron(speye(2),spones(ones(dim)));          % Specify the Jacobian pattern for implicit solvers (used for time step, not corrector step)
    end
    
    % Preparation for integration: needed for calculating the Jacobian
    Z_dim = reshape(repmat(z0_mat,dim,1),dim,n_shoot*dim);                      % This is a [dim x n_shoot*dim] matrix where each z_k is repeated dim times
    % The next matrix Delta stores the individual step widths delta_(i,j) = h*(1+abs(Z_(i,j)) for numerical differentiation
    Delta = h.*(repmat(eye(dim),1,n_shoot) + sparse(repmat(1:1:dim,1,n_shoot),1:1:n_shoot*dim,abs(s),dim,n_shoot*dim));
    Z_dim_plus  = Z_dim + Delta;                                                % This is a [dim x n_shoot*dim] matrix containing all "+ Delta" perturbed z_k
    Z_dim_minus = Z_dim - Delta;                                                % This is a [dim x n_shoot*dim] matrix containing all "- Delta" perturbed z_k
    Z_traj_plus  = zeros(dim^2,n_time,n_shoot);                                 % Stores the trajectory for the "+ delta" perturbed initial conditions
    Z_traj_minus = zeros(dim^2,n_time,n_shoot);                                 % Stores the trajectory for the "- delta" perturbed initial conditions
    Z_end_plus  = zeros(dim*dim,n_shoot);                                       % Stores the end points of the integration for the "+ delta" perturbed initial conditions
    Z_end_minus = zeros(dim*dim,n_shoot);                                       % Stores the end points of the integration for the "- delta" perturbed initial conditions
    delta_mu = h*(1+abs(mu));                                                   % Step width for numerical differentiation with respect to mu
    param_mu_plus = param;                                                      % Define a new parameter array
    param_mu_minus = param;                                                     % Define a new parameter array
    param_mu_plus{DYN.act_param}  = mu + delta_mu;                              % Update the new parameter array by mu + delta_mu
    param_mu_minus{DYN.act_param} = mu - delta_mu;                              % Update the new parameter array by mu - delta_mu
    Z_traj_mu_plus  = zeros(dim,n_time,n_shoot);                                % Stores the trajectory for "+ delta" perturbed mu
    Z_traj_mu_minus = zeros(dim,n_time,n_shoot);                                % Stores the trajectory for "- delta" perturbed mu
    Z_end_mu_plus  = zeros(dim,n_shoot);                                        % Stores the end points of the integration for the "+ delta" perturbed mu
    Z_end_mu_minus = zeros(dim,n_shoot);                                        % Stores the end points of the integration for the "- delta" perturbed mu
    Z_p_mu_plus = zeros(dim,n_time,n_shoot);                                    % Stores the trajectory at the predictor point for "+ delta" perturbed mu
    Z_p_mu_minus = zeros(dim,n_time,n_shoot);                                   % Stores the trajectory at the predictor point for "- delta" perturbed mu
    if n_auto == 0                                      % Non-autonomous case
        T_mu_plus  = 2*pi./DYN.non_auto_freq(mu+delta_mu);                      % Perturbed periodic time by "+ delta"
        T_mu_minus = 2*pi./DYN.non_auto_freq(mu-delta_mu);                      % Perturbed periodic time by "- delta"
        dT_mu_plus  = T_mu_plus/n_shoot;                                        % Perturbed time span by "+ delta" for each shooting operation (integration)
        dT_mu_minus = T_mu_minus/n_shoot;                                       % Perturbed time span by "- delta" for each shooting operation (integration)
        T_int_mu_plus  = [0:dT_mu_plus:(n_shoot-1)*dT_mu_plus;                  % Perturbed time intervals by "+ delta" for the shooting operation (integration)
                          dT_mu_plus:dT_mu_plus:n_shoot*dT_mu_plus].';    
        T_int_mu_minus = [0:dT_mu_minus:(n_shoot-1)*dT_mu_minus;                % Perturbed time intervals by "- delta" for the shooting operation (integration)
                          dT_mu_minus:dT_mu_minus:n_shoot*dT_mu_minus].';
    elseif n_auto == 1                                  % Autonomous case           
        T_int_mu_plus  = T_int;                                                 % omega cannot be the mu and thus the time intervals remain constant
        T_int_mu_minus = T_int;                                                 % omega cannot be the mu and thus the time intervals remain constant
    end

    % Integration
    % Note: The for loop is needed for non-autonomous systems since [t_start t_end] and thus Fcn(t,z0,param) is different for each time interval
    % For autonomous systems, Fcn(z,param) is independent of time t, so the integration time interval could be [0 dT(_omega)(_mu)] for all integrations
    % Using this, we could do all integrations vectorized. Instead of 3*n_shoot integrations, only 3 integrations would be necessary. However, the ...
    % vectorized integrations (only 3) can produce a small error (~ e-9 ... e-8) in Z_end. Unfortunately, this error becomes significant in the Jacobian ...
    % when approximating the derivatives and dividing by delta (~ e-8). Thus, the Jacobian can have errors ~ e-1 that eventually causes the corrector ...
    % to need more iterations. Not only does this influence the computation time negatively, it can also lead to the step control reducing the step size.
    for k=1:n_shoot
        % Create a [dim x (2*dim+2)] matrix Z0_mat storing the initial condition to integrate column-wise  
        Z0_mat = [z0_mat(:,k), Z_dim_plus(:,(k-1)*dim+1:k*dim), Z_dim_minus(:,(k-1)*dim+1:k*dim), s_p_mat(:,k)];    % Take z_k, the perturbed z_k and z_k^(P) (z_k at predictor) (need to be integrated for J)
        [~,Z] = obj.solver_function(@(t,Z) FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param)), linspace(T_int(k,1),T_int(k,2),n_time+1), Z0_mat, odeOpts_1); 
        Z_traj(:,:,k)       = Z(1:end-1,1:dim).';                   % Save the trajectory (not Z(t_end) because it will be equal to Z(t_start) of the next interval after convergence)
        Z_traj_plus(:,:,k)  = Z(1:end-1,dim+1:(dim+1)*dim).';       % Save the trajectory of the "+ delta" perturbed initial conditions
        Z_traj_minus(:,:,k) = Z(1:end-1,(dim+1)*dim+1:end-dim).';   % Save the trajectory of the "- delta" perturbed initial conditions
        Z_p(:,:,k)          = Z(1:end-1,end-dim+1:end).';           % Save the trajectory at the predictor point
        Z_end(:,k)       = Z(end,1:dim).';                          % Take the unperturbed state vectors at t_end
        Z_end_plus(:,k)  = Z(end,dim+1:(dim+1)*dim).';              % Take the "+ delta" perturbed state vectors at t_end
        Z_end_minus(:,k) = Z(end,(dim+1)*dim+1:end-dim).';          % Take the "- delta" perturbed state vectors at t_end
        % Now do a second/third integration with perturbed mu-value (needed for Jacobian - can be deactivated when it is calculated by fsolve)
        if ~(isfield(DYN.system,'first_integral') && (DYN.act_param == numel(param)))  % If the value of the first integral is NOT the continuation parameter (if it is: dg/dmu = 0 and dpc/dmu = 0)
            [~,Z_mu_plus]  = obj.solver_function(@(t,Z) FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param_mu_plus(1:end-1))), linspace(T_int_mu_plus(k,1),T_int_mu_plus(k,2),n_time+1), [z0_mat(:,k),s_p_mat(:,k)], odeOpts_2);
            Z_traj_mu_plus(:,:,k) = Z_mu_plus(1:end-1,1:dim).';         % Save the "+ delta" perturbed mu trajectory
            Z_end_mu_plus(:,k) = Z_mu_plus(end,1:dim).';                % Take the state vector of the "+ delta" perturbed mu trajectory at t_end
            Z_p_mu_plus(:,:,k) = Z_mu_plus(1:end-1,dim+1:end).';        % Save the "+ delta" perturbed mu trajectory at the predictor point
            if calc_stability                                           % Only required when using central finite difference for Jacobian
                [~,Z_mu_minus] = obj.solver_function(@(t,Z) FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param_mu_minus)), linspace(T_int_mu_minus(k,1),T_int_mu_minus(k,2),n_time+1), [z0_mat(:,k),s_p_mat(:,k)], odeOpts_2);
                Z_traj_mu_minus(:,:,k) = Z_mu_minus(1:end-1,1:dim).';   % Save the "- delta" perturbed mu trajectory
                Z_end_mu_minus(:,k) = Z_mu_minus(end,1:dim).';          % Take the state vector of the "- delta" perturbed mu trajectory at t_end
                Z_p_mu_minus(:,:,k) = Z_mu_minus(1:end-1,dim+1:end).';  % Save the "- delta" perturbed mu trajectory at the predictor point
            end
        end
    end


    % Calculate the residuum
    s_end = reshape(Z_end,dim*n_shoot,1);               % Reshape the unperturbed state vectors of Z_end (end points of shooting) to a vector
    s_perm = [s(dim+1:end); s(1:dim)];                  % Create a vector of initial conditions [z_1; z_2; ...; z_n; z_0] for residuum calculation
 
    if n_auto == 1                                      % Compute the phase condition if system is autonomous
        switch obj.phase_condition
            case 'poincare'                             
                f_s_p = reshape(Fcn(0,s_p_mat,param),dim*n_shoot,1);                                    % Evaluate the RHS at predictor point s_p
                pc = f_s_p.' * (s - reshape(s_p_mat,dim*n_shoot,1));                                    % Poincare phase condition
            case 'integral'
                f = reshape(Fcn(0,reshape(Z_traj,dim,n_shoot*n_time),param),dim*n_shoot*n_time,1);      % Evaluate the RHS at the n_shoot*n_time evaluation points
                f_p = reshape(Fcn(0,reshape(Z_p,dim,n_shoot*n_time),param),dim*n_shoot*n_time,1);       % Evaluate the RHS for the predictor point (not needed now, but later on for dpc/domega in the Jacobian)
                pc = f.' * (reshape(Z_traj,dim*n_shoot*n_time,1) - reshape(Z_p,dim*n_shoot*n_time,1));  % * T / (n_time*n_shoot); % Integral phase condition
        end
    elseif n_auto == 0                                  % Non-autonomous system
        pc = [];                                        % Set as empty in order to add "nothing" to the residuum
    end

    % If a conservative system is considered: Expand the residuum
    if isfield(DYN.system,'first_integral')
        I = DYN.system.first_integral;                              % Function of the first integral I = I(z)
        I_Z = I(z0_mat,param);                                      % Evalute the first integral for all shooting points z_i
        IC = 1/n_shoot*sum(I_Z) - param{end};                       % First Integral Constraint: I(s) = param{end} | Take the average of I_Z to get a single value for all shooting points
    else                                                            % Note: mean(I_Z) = 1/n_shoot*sum(I_Z), but the sum() function is somehow faster
        IC = [];                                                    % Non-conservative system: Set the first integral constraint IC as empty in order to add "nothing" to the residuum
    end

   res = [s_end - s_perm;  pc;  IC];                    % Assemble the residuum: Residuum of the shooting; phase condition; first integral constraint
    


    %% Jacobian J_res
    % When stability is not computed, the Jacobian is approximated using forward finite differences (sufficient)
    % When stability is computed, the Jacobian is approximated using central finite differences for higher accuracy since the Floquet multipliers are based on J
    
    % J_res = zeros(dim*n_shoot,dim*n_shoot+1);         % Needed only for testing purposes of non-autonomous systems when J is calculated by fsolve
    % J_res = zeros(dim*n_shoot+2,dim*n_shoot+2);       % Needed only for testing purposes of autonomous systems when J is calculated by fsolve
    %
    % dg/ds consists of the derivatives dz(t_(k+1),z_k,mu)/dz_k, which are [dim x dim] matrices arranged block-wise on the main diagonal of dg/ds
    % The part of dg/ds described above is named dZ_ds_mat. Furthermore, there are some -eye(dim) here and there in dg/ds. All of these are arranged in I_mat
    Z_end_dim = reshape(repmat(Z_end,dim,1),dim,n_shoot*dim);           % This is a matrix similar to Z_dim, but it stores the z(t_end,z_k)
    Z_end_plus_dim  = reshape(Z_end_plus,dim,n_shoot*dim);              % Similar to line above, but stores the ODE result of the "+ delta" perturbed z_k
    Z_end_minus_dim = reshape(Z_end_minus,dim,n_shoot*dim);             % Similar to line above, but stores the ODE result of the "- delta" perturbed z_k
    if calc_stability                                                   % Use central finite difference
        dZ_ds = (Z_end_plus_dim - Z_end_minus_dim)./repmat(2.*nonzeros(Delta)',dim,1);  % Approximate the derivatives dz(t_(k+1),z_k,mu)/dz_k
    else                                                                % Use forward finite difference
        dZ_ds = (Z_end_plus_dim - Z_end_dim)./repmat(nonzeros(Delta)',dim,1);           % Approximate the derivatives dz(t_(k+1),z_k,mu)/dz_k
    end
    dZ_ds_mat = kron(speye(n_shoot),spones(ones(dim)));                 % Create a [n_shoot*dim x n_shoot*dim] block diagonal matrix with ones(dim)
    dZ_ds_mat(logical(dZ_ds_mat)) = dZ_ds;                              % Replace all ones in dZ_ds_mat with the derivatives dz(t_(k+1),z_k,mu)/dz_k
    I_mat = spdiags(ones(n_shoot*dim,1),dim,n_shoot*dim,n_shoot*dim);   % Create a matrix where "1" are placed on the "dim"-th upper right secondary diagonal
    I_mat(end-dim+1:end,1:dim) = eye(dim);                              % eye(dim) needs to be added in the bottom left corner
    dg_ds = dZ_ds_mat - I_mat;                                          % dg/ds

    
    % Now calculate dg/dmu
    if isfield(DYN.system,'first_integral') && (DYN.act_param == numel(param))  % If the value of the first integral is the continuation parameter
        dg_dmu = zeros(dim*n_shoot,1);                                  % dg/dmu is zero since g is not explicitly dependent on I = mu
    else
        if calc_stability                                               % Use central finite difference
            dg_dmu = reshape( (Z_end_mu_plus - Z_end_mu_minus) ./ (2*delta_mu) , dim*n_shoot, 1);
        else                                                            % Use forward finite difference
            dg_dmu = reshape( (Z_end_mu_plus - Z_end) ./ delta_mu , dim*n_shoot, 1);
        end
    end

   
    % Autonomous system: Calculate the derivative of g with respect to omega and the derivative of the phase condition
    if n_auto == 1

        % dg/domega:
        f_end = Fcn(0,Z_end,param);                                                     % Evaluate RHS at end points of integration
        dg_domega = - 2*pi / (n_shoot * omega^2) .* reshape(f_end,dim*n_shoot,1);       % dg/domega

        % dpc/dy:
        switch obj.phase_condition
            
            case 'poincare'                     % Poincare phase condition
                dpc_ds = f_s_p.';                                                                               % dpc/ds
                dpc_domega = 0;                                                                                 % dpc/domega = 0 (phase condition is independent of the frequency)
                f_mu_plus = reshape(Fcn(0,s_p_mat,param_mu_plus),dim*n_shoot,1);                                % Evaluate the RHS at predictor point s_p with perturbed mu
                if calc_stability                                                                               % Use central finite difference
                    f_mu_minus = reshape(Fcn(0,s_p_mat,param_mu_minus),dim*n_shoot,1);                          % Evaluate the RHS at predictor point s_p with perturbed mu
                    dpc_dmu = (f_mu_plus - f_mu_minus).' * (s - reshape(s_p_mat,dim*n_shoot,1)) / (2*delta_mu); % dpc/dmu
                else                                                                                            % Use forward finite difference
                    dpc_dmu = (f_mu_plus - f_s_p).' * (s - reshape(s_p_mat,dim*n_shoot,1)) / delta_mu;          % dpc/dmu
                end
            
            case 'integral'                     % Integral phase condition
                dpc_ds = zeros(1,dim*n_shoot);                                      % Initialise
                dpc_domega = 0;                                                     % Initialise -> dpc_domega is iteratively computed and added up in a loop
                for k = 1:n_shoot                                                   % Loop over the shooting points / interval -> Implementation without very difficult and complicated
                    Delta_k = repmat(Delta(:,(k-1)*dim+1:k*dim),1,n_time);          % Part of matrix Delta that belongs to k-th shooting point and repeated n_time times -> size: [dim x dim*n_time]
                    f_k = f((k-1)*dim*n_time+1:k*dim*n_time);                       % Part of f (Fcn evaluated for all ~100 evaluation points) that belongs to k-th interval -> size: [dim*n_time x 1]
                    f_p_k = f_p((k-1)*dim*n_time+1:k*dim*n_time);                   % Like f above, but for the predictor point
                    Z_traj_k = reshape(Z_traj(:,:,k),dim*n_time,1);                 % Solution trajectory within the k-th interval -> size: [dim*n_time x 1]
                    Z_p_k = reshape(Z_p(:,:,k),dim*n_time,1);                  	    % Predictor trajectory within the k-th interval -> size: [dim*n_time x 1]
                    % We need to compute some derivatives w.r.t. the evaluation points z_j = z(t_j,z_k), where the index j describes the evaluation points within the k-th interval -> j = 1, ..., n_time and z_{j=1} = z_k
                    % Therefore: The z_j needs to be perturbed like the shooting points above. We also need a corresponding Delta matrix storing the individual step widths delta_(j,m) = h*(1+abs(Z_(j,m)) for num. diff.
                    Z_dim_j = reshape(repmat(Z_traj(:,:,k),dim,1),dim,n_time*dim);  % This is a [dim x n_time*dim] matrix where each z_j is repeated dim times
                    Delta_j = h.*(repmat(eye(dim),1,n_time) + sparse(repmat(1:1:dim,1,n_time),1:1:n_time*dim,abs(Z_traj_k),dim,n_time*dim));
                    Z_dim_j_plus  = Z_dim_j + Delta_j;                              % This is a [dim x n_time*dim] matrix containing all "+ Delta" perturbed z_j
                    Z_dim_j_minus = Z_dim_j - Delta_j;                              % This is a [dim x n_time*dim] matrix containing all "- Delta" perturbed z_j
                    if calc_stability
                        % Approximate the derivatives df(z_j(z_k),mu)/dz_k. The matrices are arranged as: [df(z_0(z_k),mu)/dz_k, df(z_1(z_k),mu)/dz_k, ...] = [dim x dim*n_time]
                        df_dz_k = (Fcn(0,reshape(Z_traj_plus(:,:,k),dim,dim*n_time),param) - Fcn(0,reshape(Z_traj_minus(:,:,k),dim,dim*n_time),param)) ./ repmat(2.*nonzeros(Delta_k)',dim,1);
                        % Approximate the derivatives df(z_j,mu)/dz_j (matrices arranged like above)
                        df_dz_j = (Fcn(0,Z_dim_j_plus,param) - Fcn(0,Z_dim_j_minus,param)) ./ repmat(2.*nonzeros(Delta_j)',dim,1);
                        % Approximate the derivatives dz_j(z_k)/dz_k (matrices arranged like above)
                        dz_j_dz_k = (reshape(Z_traj_plus(:,:,k),dim,dim*n_time) - reshape(Z_traj_minus(:,:,k),dim,dim*n_time)) ./ repmat(2.*nonzeros(Delta_k)',dim,1);
                    else
                        % Approximate the derivatives df(z_j(z_k),mu)/dz_k. The matrices are arranged as: [df(z_0(z_k),mu)/dz_k, df(z_1(z_k),mu)/dz_k, ...] = [dim x dim*n_time]
                        df_dz_k = (Fcn(0,reshape(Z_traj_plus(:,:,k),dim,dim*n_time),param) - Fcn(0,Z_dim_j,param)) ./ repmat(nonzeros(Delta_k)',dim,1);
                        % Approximate the derivatives df(z_j,mu)/dz_j (matrices arranged like above)
                        df_dz_j = (Fcn(0,Z_dim_j_plus,param) - Fcn(0,Z_dim_j,param)) ./ repmat(nonzeros(Delta_j)',dim,1);
                        % Approximate the derivatives dz_j(z_k)/dz_k (matrices arranged like above)
                        dz_j_dz_k = (reshape(Z_traj_plus(:,:,k),dim,dim*n_time) - reshape(Z_dim_j,dim,dim*n_time)) ./ repmat(nonzeros(Delta_k)',dim,1);
                    end
                    df_dz_k = reshape(permute(reshape(df_dz_k,dim,dim,n_time),[1 3 2]),dim*n_time,dim);         % Reshape matrices [M1, M2, ... , Md] ([dim x dim*n_time]) into [M1; M2; ... ; Md] ([dim*n_time x dim])
                    dz_j_dz_k = reshape(permute(reshape(dz_j_dz_k,dim,dim,n_time),[1 3 2]),dim*n_time,dim);     % Reason: We need this particular structure to compute dpc/dz_k by matrix multiplication without a loop
                    % dpc_ds:
                    % Do a loop over the n_time evaluation points (t_j,z_j) -> only necessary to compute df_dz_k = df_dz_j * dz_j_dz_k or to compute dpc_ds via direct numerical differentiation
                    % for j = 1:n_time
                        % df_dz_k((j-1)*dim+1:j*dim,:) = df_dz_j((j-1)*dim+1:j*dim,:) * dz_j_dz_k((j-1)*dim+1:j*dim,:);       % Alternative way of computing df(z_j(z_k),mu)/dz_k -> quasi identical to method above
                        % Alternative method to compute dpc_ds: "Direct" numerical differentiation of the phase condition without exploiting product (and chain) rule -> slower than other method
                        % Z_k_j_plus = reshape(Z_traj_plus(:,j,k),dim,dim);         % [dim x dim] matrix of the d different z_j after perturbing z_k by "+ delta" in each dimension
                        % Z_k_j_minus = reshape(Z_traj_minus(:,j,k),dim,dim);       % [dim x dim] matrix of the d different z_j after perturbing z_k by "- delta" in each dimension
                        % Z_j_mat = Z_dim_j(:,(j-1)*dim+1:j*dim);                   % [dim x dim] matrix where the unperturbed z_j is repeated d times
                        % Z_p_j = repmat(Z_p_traj(:,j,k),1,dim);                    % [dim x dim] matrix where z_p(t_j) (z_j at predictor point) is repeated d times
                        % if calc_stability                                         % If this method is used: temp = zeros(1,dim) needs to be initialised before the loop is executed
                        %     dpc_dz_k = dpc_dz_k + diag(Fcn(0,Z_k_j_plus,param).' * (Z_k_j_plus - Z_p_j) - Fcn(0,Z_k_j_minus,param).' * (Z_k_j_minus - Z_p_j)).' ./ (2.*nonzeros(Delta_k(:,1:dim))).';
                        % else
                        %     dpc_dz_k = dpc_dz_k + diag(Fcn(0,Z_k_j_plus,param).' * (Z_k_j_plus - Z_p_j) - Fcn(0,Z_j_mat,param).' * (Z_j_mat - Z_p_j)).' ./ nonzeros(Delta_k(:,1:dim)).';
                        % end
                    % end
                    % dpc_ds(1,(k-1)*dim+1:k*dim) = dpc_dz_k;                                                   % dpc/dz_k when using the "direct" numerical differentiation
                    dpc_ds(1,(k-1)*dim+1:k*dim) = (Z_traj_k - Z_p_k).' * df_dz_k + f_k.' * dz_j_dz_k;           % dpc/dz_k
                    % dpc_domega:
                    df_dz_j_mat = kron(speye(n_time),spones(ones(dim)));            % Create a [n_time*dim x n_time*dim] block diagonal matrix with ones(dim)
                    df_dz_j_mat(logical(df_dz_j_mat)) = df_dz_j;                    % Replace all ones in df_dz_j_mat with the derivatives df(z_j,mu)/dz_j
                    j_vec = reshape(repmat(0:n_time-1,dim,1),dim*n_time,1);         % Create a vector containing dim times 0, dim times 1, ..., dim times (n_time-1)
                    dpc_domega = dpc_domega - 2*pi/(n_time*n_shoot*omega^2) * ((df_dz_j_mat*f_k).' * (j_vec.*(Z_traj_k - Z_p_k)) + f_k.'*(j_vec.*(f_k-f_p_k)));  
                end
                % dpc_dmu:
                if isfield(DYN.system,'first_integral') && (DYN.act_param == numel(param))  % If the value of the first integral is the continuation parameter
                    dpc_dmu = 0;                                                    % dpc/dmu is zero since pc is not explicitly dependent on I = mu
                else
                    if calc_stability
                        df_dmu = (Fcn(0,reshape(Z_traj_mu_plus,dim,n_time*n_shoot),param_mu_plus) - Fcn(0,reshape(Z_traj_mu_minus,dim,n_time*n_shoot),param_mu_minus)) ./ (2*delta_mu);
                        dz_j_dmu = (reshape(Z_traj_mu_plus,dim,n_time*n_shoot) - reshape(Z_traj_mu_minus,dim,n_time*n_shoot)) ./ (2*delta_mu);
                        dz_p_dmu = (reshape(Z_p_mu_plus,dim,n_time*n_shoot) - reshape(Z_p_mu_minus,dim,n_time*n_shoot)) ./ (2*delta_mu);
                    else
                        df_dmu = (Fcn(0,reshape(Z_traj_mu_plus,dim,n_time*n_shoot),param_mu_plus) - Fcn(0,reshape(Z_traj,dim,n_time*n_shoot),param)) ./ delta_mu;
                        dz_j_dmu = (reshape(Z_traj_mu_plus,dim,n_time*n_shoot) - reshape(Z_traj,dim,n_time*n_shoot)) ./ delta_mu;
                        dz_p_dmu = (reshape(Z_p_mu_plus,dim,n_time*n_shoot) - reshape(Z_p,dim,n_time*n_shoot)) ./ delta_mu;
                    end
                    dpc_dmu = reshape(df_dmu,dim*n_time*n_shoot,1).' * reshape(Z_traj - Z_p,dim*n_time*n_shoot,1) + f.' * reshape(dz_j_dmu - dz_p_dmu,dim*n_time*n_shoot,1);
                end
        end

    % Non-Autonomous system: Set dg/domega and dpc/dy as empty in order to add "nothing" to the Jacobian matrix below
    elseif n_auto == 0
        dg_domega = [];  dpc_ds = [];  dpc_domega = [];  dpc_dmu = [];      
    end
    

    % Conservative system: Calculate dIC/dy -> GEHT NOCH NICHT, DRÜBER SCHAUEN (VARIABLENNAMEN UND SO)
    if isfield(DYN.system,'first_integral')
        if calc_stability                                                                                           % Use central finite difference
            dIC_ds = 1/n_shoot .* (I(Z_dim_plus,param) - I(Z_dim_minus,param)) ./ (2.*nonzeros(Delta)');
        else                                                                                                        % Use forward finite difference
            dIC_ds = 1/n_shoot .* (I(Z_dim_plus,param) - reshape(repmat(I_Z,dim,1),1,n_shoot*dim)) ./ nonzeros(Delta)';
        end
        if n_auto == 1;     dIC_domega = 0;                                                                         % Autonomous system: dIC/domega is also required
        else;               dIC_domega = [];                                                                        % Non-Autonomous system: Set dIC/domega as empty in order to add "nothing" to the Jacobian matrix below
        end
        if DYN.act_param == numel(param)                                                                            % If param{end} = first integral is the continuation parameter
            dIC_dmu = - 1;                                                                                          % I(Z,param) is independent of mu and dparam{end}/dmu = 1
        else                                                                                                        % A parameter of the RHS is the continuation parameter: In this case, I can be dependent on mu
            if calc_stability                                                                                       % Use central finite difference
                dIC_dmu = 1/n_shoot*sum(I(z0_mat,param_mu_plus) - I(z0_mat,param_mu_minus)) / (2*delta_mu);
            else                                                                                                    % Use forward finite difference
                dIC_dmu = 1/n_shoot*sum(I(z0_mat,param_mu_plus) - I_Z) / delta_mu;
            end
        end
    % Non-Conservative system: Set dIC/dy as empty in order to add "nothing" to the Jacobian matrix below
    else
        dIC_ds = [];  dIC_domega = [];  dIC_dmu = [];
    end

    % Build the Jacobian matrix J_res
    J_res = [dg_ds,  dg_domega,  dg_dmu;                                % Build the Jacobian matrix
             dpc_ds, dpc_domega, dpc_dmu;                               % Build the Jacobian matrix
             dIC_ds, dIC_domega, dIC_dmu];                              % Build the Jacobian matrix
    %}


end



%% A function wrapper is needed for the ODE-solver to be able to integrate multiple initial conditions simultaneously
%
% @t:   time
% @Z:   hyper-vector of state vectors
% @dim: Dimension of the system
% @Fcn: RHS of ODE

function dZdt = FCN_wrapper(t,Z,dim,Fcn)

    n = numel(Z)/dim;                       % Number of state vectors
    z_mat = reshape(Z,dim,n);               % Reshape the hyper-vector to a matrix where the initial conditions are arranged column by column
    dzdt = Fcn(t,z_mat);                    % Evaluate the RHS vectorized
    dZdt = reshape(dzdt,dim*n,1);           % Reshape the [dim x n] - matrix dzdt to a hyper-vector

end