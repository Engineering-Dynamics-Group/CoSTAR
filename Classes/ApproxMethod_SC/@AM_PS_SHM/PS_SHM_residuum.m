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
    n_time = obj.n_time;                                % Number of time evaluation points in each shooting interval for the integral phase condition
    s = y(1:(end-1-n_auto));                            % Method solution vector (hyper-vector of shooting points [z_0; z_1; ...; z_n])
    s0 = obj.iv(1:end-n_auto);                          % Method solution vector at the predictor point (obj.iv(end) is the autonomous frequency if present)
    mu = y(end);                                        % Continuation parameter
    Z0 = obj.Z0;                                        % Get the trajectory of the preceding solution for evaluating the integral phase condition
    
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
    Z_traj = zeros(dim,n_time,n_shoot);                 % Stores the trajectory for evaluating the integral phase condition
    Z_end = zeros(dim,n_shoot);                         % Stores the end points of integration (could be stored in Z_traj as well, but using Z_end is more convenient)
    odeOpts_1 = obj.odeOpts;                            % Use these options for the integrations with perturbed z_i (because that is where the FCN_wrapper is used)
    if strcmpi(obj.solver,'ode15s') || strcmpi(obj.solver,'ode23s') || strcmpi(obj.solver,'ode23t') || strcmpi(obj.solver,'ode23tb')
        odeOpts_1.JPattern = kron(speye(2*dim+1),spones(ones(dim)));    % Specify the Jacobian pattern for implicit solvers (used for time step, not corrector step)
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
        Z0_mat = [z0_mat(:,k), Z_dim_plus(:,(k-1)*dim+1:k*dim), Z_dim_minus(:,(k-1)*dim+1:k*dim)];      % Take z_k and the perturbed z_k (need to be integrated for J)
        [~,Z] = obj.solver_function(@(t,Z) FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param)), linspace(T_int(k,1),T_int(k,2),n_time+1), Z0_mat, odeOpts_1); 
        Z_traj(:,:,k)       = Z(1:end-1,1:dim).';                   % Save the trajectory (not Z(t_end) because it will be equal to Z(t_start) of the next interval after convergence)
        Z_traj_plus(:,:,k)  = Z(1:end-1,dim+1:(dim+1)*dim).';       % Save the trajectory of the "+ delta" perturbed initial conditions
        Z_traj_minus(:,:,k) = Z(1:end-1,end-dim*dim+1:end).';       % Save the trajectory of the "- delta" perturbed initial conditions
        Z_end(:,k)       = Z(end,1:dim).';                          % Take the unperturbed state vectors at t_end
        Z_end_plus(:,k)  = Z(end,dim+1:(dim+1)*dim).';              % Take the "+ delta" perturbed state vectors at t_end
        Z_end_minus(:,k) = Z(end,end-dim*dim+1:end).';              % Take the "- delta" perturbed state vectors at t_end
        % Now do a second/third integration with perturbed mu-value (needed for Jacobian - can be deactivated when it is calculated by fsolve)
        if ~(isfield(DYN.system,'first_integral') && (DYN.act_param == numel(param)))  % If the value of the first integral is NOT the continuation parameter (if it is: dg/dmu = 0 and dpc/dmu = 0)
            [~,Z_mu_plus]  = obj.solver_function(@(t,Z) FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param_mu_plus)), linspace(T_int_mu_plus(k,1),T_int_mu_plus(k,2),n_time+1), z0_mat(:,k), obj.odeOpts);
            Z_traj_mu_plus(:,:,k) = Z_mu_plus(1:end-1,1:dim).';         % Save the "+ delta" perturbed mu trajectory
            Z_end_mu_plus(:,k) = Z_mu_plus(end,1:dim).';                % Take the state vector of the "+ delta" perturbed mu trajectory at t_end
            if calc_stability                                           % Only required when using central finite difference for Jacobian
                [~,Z_mu_minus] = obj.solver_function(@(t,Z) FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param_mu_minus)), linspace(T_int_mu_minus(k,1),T_int_mu_minus(k,2),n_time+1), z0_mat(:,k), obj.odeOpts);
                Z_traj_mu_minus(:,:,k) = Z_mu_minus(1:end-1,1:dim).';   % Save the "- delta" perturbed mu trajectory
                Z_end_mu_minus(:,k) = Z_mu_minus(end,1:dim).';          % Take the state vector of the "- delta" perturbed mu trajectory at t_end
            end
        end
    end


    % Calculate the residuum
    s_end = reshape(Z_end,dim*n_shoot,1);               % Reshape the unperturbed state vectors of Z_end (end points of shooting) to a vector
    s_perm = [s(dim+1:end); s(1:dim)];                  % Create a vector of initial conditions [z_1; z_2; ...; z_n; z_0] for residuum calculation
 
    if n_auto == 1                                                                                      % Compute the phase condition if system is autonomous
        switch obj.phase_condition
            case 'poincare'
                z_0 = s(1:dim);                                                                         % State space vector at theta = 0
                z0_0 = s0(1:dim);                                                                       % State space vector at the preceding point for theta = 0
                f_z0_0 = Fcn(0,z0_0,param);                                                             % Evaluate Fcn for the phase condition
                pc = f_z0_0' * (z_0 - z0_0);                                                            % Poincare phase condition
            case 'integral'
                obj.Z_traj = Z_traj;                                                                    % Save the trajectory so that it can be used directly as Z0 when the next solution is computed (see IF_up_res_data)
                f = reshape(Fcn(0,reshape(Z_traj,dim,n_shoot*n_time),param),dim*n_shoot*n_time,1);      % Evaluate the RHS at the n_shoot*n_time evaluation points
                f0 = reshape(Fcn(0,reshape(Z0,dim,n_shoot*n_time),param),dim*n_shoot*n_time,1);         % Evaluate the RHS for the predictor point
                pc = f0.' * (reshape(Z_traj,dim*n_shoot*n_time,1) - reshape(Z0,dim*n_shoot*n_time,1));  % * T / (n_time*n_shoot); % Integral phase condition
        end
    else                                                % Non-autonomous system OR autonomous conservative system without phase condition
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

   
    % Autonomous system: Calculate dg/domega and dpc/dy
    if n_auto == 1

        % dg/domega:
        f_end = Fcn(0,Z_end,param);                                                     % Evaluate RHS at end points of integration
        dg_domega = - 2*pi / (n_shoot * omega^2) .* reshape(f_end,dim*n_shoot,1);       % dg/domega

        % dpc/dy:
        switch obj.phase_condition 
            case 'poincare'                     % Poincare phase condition
                dpc_ds = [f_z0_0', zeros(1,(n_shoot-1)*dim)];                                                           % dpc/ds of Poincare phase condition
                dpc_domega = 0;                                                                                         % dpc/domega = 0 (phase condition is independent of the frequency)
                if calc_stability                                                                                       % Use central finite difference
                    dpc_dmu = (Fcn(0,z0_0,param_mu_plus) - Fcn(0,z0_0,param_mu_minus)).' * (z_0 - z0_0) / (2*delta_mu); % dpc/dmu
                else                                                                                                    % Use forward finite difference
                    dpc_dmu = (Fcn(0,z0_0,param_mu_plus) - f_z0_0).' * (z_0 - z0_0) / delta_mu;                         % dpc/dmu
                end
            case 'integral'                     % Integral phase condition
                % Calculating dpc_ds and dpc_domega with a loop: The loop is not required, but this code is easier to understand than the vectorized variant below
                % dpc_ds = zeros(1,dim*n_shoot);                                      % Initialise
                % dpc_domega = 0;                                                     % Initialise -> dpc_domega is iteratively computed and added up in a loop
                % for k = 1:n_shoot
                %     Delta_k = repmat(Delta(:,(k-1)*dim+1:k*dim),1,n_time);          % Part of matrix Delta that belongs to k-th shooting point and repeated n_time times -> size: [dim x dim*n_time]
                %     f_k = f((k-1)*dim*n_time+1:k*dim*n_time);                       % Part of f (Fcn evaluated for all ~100 evaluation points) that belongs to k-th interval -> size: [dim*n_time x 1]
                %     f0_k = f0((k-1)*dim*n_time+1:k*dim*n_time);                     % Part of f (Fcn evaluated for all ~100 evaluation points) that belongs to k-th interval -> size: [dim*n_time x 1]
                %     Z_dim_j = repmat(Z_traj(:,:,k),dim,1);                          % This is a [dim*dim x n_time] matrix where each z_j is repeated dim times below
                %     if calc_stability                                               % Approximate the derivatives dz_j(z_k)/dz_k (matrices arranged like above)
                %         dz_j_dz_k = reshape(Z_traj_plus(:,:,k) - Z_traj_minus(:,:,k),dim,dim*n_time) ./ repmat(2.*nonzeros(Delta_k)',dim,1);
                %     else
                %         dz_j_dz_k = reshape(Z_traj_plus(:,:,k) - Z_dim_j,dim,dim*n_time) ./ repmat(nonzeros(Delta_k)',dim,1);
                %     end
                %     dz_j_dz_k = reshape( permute( reshape(dz_j_dz_k,dim,dim,n_time), [1 3 2] ), dim*n_time, dim);   % We need this particular structure to compute dpc/dz_k by matrix multiplication without a loop
                %     % dz_j_dz_k now stores the n_time [dim x dim] matrices dz_j/dz_k below each other (before the line above, these matrices were stored right next to each other)
                %     dpc_ds(1,(k-1)*dim+1:k*dim) = f0_k.' * dz_j_dz_k;                                 % dpc/dz_k
                %     dpc_domega = dpc_domega - 2*pi/(n_time*n_shoot*omega^2) * f0_k.'*(j_vec.*f_k);    % dpc/domega
                % end
                % Calculating dpc_ds and dpc_domega without a loop:
                Delta_exp = repmat(reshape(Delta,dim,dim,n_shoot),1,n_time,1);      % We need the Delta matrix (n_time)-times to compute dz_j_dz_k: Delta = [dim x dim*n_shoot] -> Delta_exp = [dim x dim*n_time x n_shoot]
                Z_traj_dim = repmat(Z_traj,dim,1,1);                                % This is needed to compute dz_j_dz_k via forward difference
                if calc_stability                                                   % Approximate the derivatives dz_j(z_k)/dz_k
                    dz_j_dz_k = reshape(Z_traj_plus - Z_traj_minus,dim,dim*n_time,n_shoot) ./ repmat(2.*reshape(nonzeros(Delta_exp),1,dim*n_time,n_shoot),dim,1);
                else
                    dz_j_dz_k = reshape(Z_traj_plus - Z_traj_dim,dim,dim*n_time,n_shoot) ./ repmat(reshape(nonzeros(Delta_exp),1,dim*n_time,n_shoot),dim,1);
                end
                dz_j_dz_k = reshape( permute( reshape(dz_j_dz_k,dim,dim,n_time,n_shoot), [1 3 2 4] ), dim*n_time, dim*n_shoot );    % We need this particular structure for the multiplikation with f0
                % dz_j_dz_k stores the n_time*n_shoot [dim x dim] matrices dz_j/dz_k in the k-th "column" and the j-th "row". The reshape(f0,...).' below stores in the k-th row the f(z0_{k,j}) for all j 
                f0_dz_j_dz_k = reshape(f0,dim*n_time,n_shoot).' * dz_j_dz_k;      % This product evaluates all f(z0_{k,j}).' * dz_j/dz_k (for all combinations of the two k). However, we only need products of matching k
                eval_mat = logical( kron(speye(n_shoot),ones(1,dim)) );             % Evaluation matrix to obtain the products of matching k. They can be found on the "pseudo-diagonal" defined by the kron()-product
                dpc_ds = reshape(f0_dz_j_dz_k(eval_mat),1,dim*n_shoot);            % dpc/ds | The reshape is needed because a row vector is returned for n_shoot = 1, while a column vector is returned for > 1
                j_vec = reshape(repmat(0:n_time-1,dim,1),dim*n_time,1);             % Create a vector containing dim times 0, dim times 1, ..., dim times (n_time-1)
                J_vec = repmat(j_vec,n_shoot,1);                                    % j_vec is repeated n_shoot times since j_vec is only for the k-th shooting point (but it is constant for all k)
                dpc_domega = - 2*pi/(n_time*n_shoot*omega^2) * f0.'*(J_vec.*f);    % dpc_domega
                % dpc_dmu:
                if isfield(DYN.system,'first_integral') && (DYN.act_param == numel(param))  % If the value of the first integral is the continuation parameter
                    dpc_dmu = 0;                                                    % dpc/dmu is zero since pc is not explicitly dependent on I = mu
                else
                    if calc_stability
                        df0_dmu = (Fcn(0,reshape(Z0,dim,n_time*n_shoot),param_mu_plus) - Fcn(0,reshape(Z0,dim,n_time*n_shoot),param_mu_minus)) ./ (2*delta_mu);
                        dz_j_dmu = (reshape(Z_traj_mu_plus,dim,n_time*n_shoot) - reshape(Z_traj_mu_minus,dim,n_time*n_shoot)) ./ (2*delta_mu);
                    else
                        df0_dmu = (Fcn(0,reshape(Z0,dim,n_time*n_shoot),param_mu_plus) - Fcn(0,reshape(Z0,dim,n_time*n_shoot),param)) ./ delta_mu;
                        dz_j_dmu = (reshape(Z_traj_mu_plus,dim,n_time*n_shoot) - reshape(Z_traj,dim,n_time*n_shoot)) ./ delta_mu;
                    end
                    dpc_dmu = reshape(df0_dmu,dim*n_time*n_shoot,1).' * reshape(Z_traj - Z0,dim*n_time*n_shoot,1) + f0.' * reshape(dz_j_dmu,dim*n_time*n_shoot,1);
                end
        end
    % Phase condition is missing: Set dg/domega and dpc/dy as empty in order to add "nothing" to the Jacobian matrix below
    else
        dg_domega = [];  dpc_ds = [];  dpc_domega = [];  dpc_dmu = [];      
    end
    

    % Conservative system: Calculate dIC/dy
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