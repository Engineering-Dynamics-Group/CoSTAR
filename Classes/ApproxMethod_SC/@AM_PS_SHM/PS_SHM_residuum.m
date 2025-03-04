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
    s = y(1:(end-1-n_auto));                            % Method solution vector (hyper-vector of shooting points [z_0; z_1; ...; z_n])
    mu = y(end);                                        % Continuation parameter

    if n_auto == 0
        omega = DYN.non_auto_freq(mu);                  % Angular frequency (non-autonomous system)
    elseif n_auto == 1
        omega = y(end-1);                               % Angular frequency (autonomous system)
        s_p = obj.iv(1:end-1);                          % Solution vector s at the predictor point (obj.iv(end) is the autonomous frequency at the predictor point)
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
    Z_end = zeros(dim,n_shoot);                         % Initialize array for end points of shooting
    Z_end_plus  = zeros(dim*dim,n_shoot);               % Stores the end points of shooting with "+ delta" perturbed initial conditions
    Z_end_minus = zeros(dim*dim,n_shoot);               % Stores the end points of shooting with "- delta" perturbed initial conditions

    % Preparation for integration: needed for calculating the Jacobian
    Z_dim = reshape(repmat(z0_mat,dim,1),dim,n_shoot*dim);                      % This is a [dim x n_shoot*dim] matrix where each z_i is repeated dim times
    % The next matrix Delta stores the individual step widths delta_(i,j) = h*(1+abs(Z_(i,j)) for numerical differentiation
    Delta = h.*(repmat(eye(dim),1,n_shoot) + sparse(repmat(1:1:dim,1,n_shoot),1:1:n_shoot*dim,abs(s),dim,n_shoot*dim));
    Z_dim_plus  = Z_dim + Delta;                                                % This is a [dim x n_shoot*dim] matrix containing all "+ Delta" perturbed z_i
    Z_dim_minus = Z_dim - Delta;                                                % This is a [dim x n_shoot*dim] matrix containing all "- Delta" perturbed z_i
    delta_mu = h*(1+abs(mu));                                                   % Step width for numerical differentiation with respect to mu
    param_mu_plus = param;                                                      % Define a new parameter array
    param_mu_minus = param;                                                     % Define a new parameter array
    param_mu_plus{DYN.act_param}  = mu + delta_mu;                              % Update the new parameter array by mu + delta_mu
    param_mu_minus{DYN.act_param} = mu - delta_mu;                              % Update the new parameter array by mu - delta_mu
    Z_end_mu_plus  = zeros(dim,n_shoot);                                        % Stores the end points of shooting with "+ delta" perturbed mu
    Z_end_mu_minus = zeros(dim,n_shoot);                                        % Stores the end points of shooting with "- delta" perturbed mu
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
        % Create a [dim x (2*dim+1)] matrix Z0_mat storing the initial condition to integrate column-wise
        Z0_mat = [z0_mat(:,k), Z_dim_plus(:,(k-1)*dim+1:k*dim), Z_dim_minus(:,(k-1)*dim+1:k*dim)];   % Take z_i and the perturbed z_i (need to be integrated for J)
        [~,Z] = obj.solver_function(@(t,Z)FCN_wrapper(t,Z,dim,@(t,z)Fcn(t,z,param)), T_int(k,:), Z0_mat, obj.odeOpts); 
        Z_end(:,k) = Z(end,1:dim).';                                % Take the unperturbed state vectors at t_end
        Z_end_plus(:,k)  = Z(end,dim+1:(dim+1)*dim).';              % Take the "+ delta" perturbed state vectors at t_end
        Z_end_minus(:,k) = Z(end,(dim+1)*dim+1:end).';              % Take the "- delta" perturbed state vectors at t_end
        % Now do a second/third integration with perturbed mu-value (needed for Jacobian - can be deactivated when it is calculated by fsolve)
        [~,Z_mu_plus]  = obj.solver_function(@(t,z) Fcn(t,z,param_mu_plus), T_int_mu_plus(k,:), z0_mat(:,k), obj.odeOpts); 
        Z_end_mu_plus(:,k) = Z_mu_plus(end,:).';                    % Take the "+ delta" perturbed mu state vectors at t_end
        if calc_stability                                           % Only required when using central finite difference for Jacobian
            [~,Z_mu_minus] = obj.solver_function(@(t,z) Fcn(t,z,param_mu_minus), T_int_mu_minus(k,:), z0_mat(:,k), obj.odeOpts); 
            Z_end_mu_minus(:,k) = Z_mu_minus(end,:).';              % Take the "- delta" perturbed mu state vectors at t_end
        end
    end


    % Calculate the residuum
    s_end = reshape(Z_end,dim*n_shoot,1);               % Reshape the unperturbed state vectors of Z_end (end points of shooting) to a vector
    s_perm = [s(dim+1:end); s(1:dim)];                  % Create a vector of initial conditions [z_1; z_2; ...; z_n; z_0] for residuum calculation
    if n_auto == 0
        res = s_end - s_perm;                           % Residuum
    elseif n_auto == 1
        f_p = reshape(Fcn(0,reshape(s_p,dim,n_shoot),param),dim*n_shoot,1);     % Evaluate the RHS at predictor point s_p
        res = [s_end - s_perm;                          % Residuum
               f_p.' * (s - s_p)];                      % Poincare phase condition ( = pc)
    end
    


    %% Jacobian J_res
    % When stability is not computed, the Jacobian is approximated using forward finite differences (sufficient)
    % When stability is computed, the Jacobian is approximated using central finite differences for higher accuracy since the Floquet multipliers are based on J
    
    % J_res = zeros(dim*n_shoot,dim*n_shoot+1);         % Needed only for testing purposes of non-autonomous systems when J is calculated by fsolve
    % J_res = zeros(dim*n_shoot+1,dim*n_shoot+2);       % Needed only for testing purposes of autonomous systems when J is calculated by fsolve
    %
    % dg/ds consists of the derivatives dz(t_(i+1),z_i,mu)/dz_i, which are [dim x dim] matrices arranged block-wise on the main diagonal of dg/ds
    % The part of dg/ds described above is named dZ_ds_mat. Furthermore, there are some -eye(dim) here and there in dg/ds. All of these are arranged in I_mat
    Z_end_dim = reshape(repmat(Z_end,dim,1),dim,n_shoot*dim);           % This is a matrix similar to Z_dim, but it stores the z(t_end,z_i)
    Z_end_plus_dim  = reshape(Z_end_plus,dim,n_shoot*dim);              % Similar to line above, but stores the ODE result of the "+ delta" perturbed z_i
    Z_end_minus_dim = reshape(Z_end_minus,dim,n_shoot*dim);             % Similar to line above, but stores the ODE result of the "- delta" perturbed z_i
    if calc_stability                                                   % Use central finite difference
        dZ_ds = (Z_end_plus_dim - Z_end_minus_dim)./repmat(2.*nonzeros(Delta)',dim,1);  % Approximate the derivatives dz(t_(i+1),z_i,mu)/dz_i
    else                                                                % Use forward finite difference
        dZ_ds = (Z_end_plus_dim - Z_end_dim)./repmat(nonzeros(Delta)',dim,1);           % Approximate the derivatives dz(t_(i+1),z_i,mu)/dz_i
    end
    dZ_ds_mat = kron(speye(n_shoot),spones(ones(dim)));                 % Create a [n_shoot*dim x n_shoot*dim] block diagonal matrix with ones(dim)
    dZ_ds_mat(logical(dZ_ds_mat)) = dZ_ds;                              % Replace all ones in dZ_ds_mat with the derivatives dz(t_(i+1),z_i,mu)/dz_i
    I_mat = kron(spdiags([0 1],0:1,n_shoot,n_shoot),eye(dim));          % Create a matrix where eye(dim) is placed on the first upper right secondary diagonal
    I_mat(end-dim+1:end,1:dim) = eye(dim);                              % eye(dim) needs to be added in the bottom left corner
    dg_ds = dZ_ds_mat - I_mat;                                          % dg/ds
    
    
    % Now calculate dg/dmu
    if calc_stability                                                   % Use central finite difference
        dg_dmu = reshape( (Z_end_mu_plus - Z_end_mu_minus) ./ (2*delta_mu) , dim*n_shoot, 1);
    else                                                                % Use forward finite difference
        dg_dmu = reshape( (Z_end_mu_plus - Z_end) ./ delta_mu , dim*n_shoot, 1);
    end


    % Build the Jacobian matrix J_res
    if n_auto == 0                                      % Non-autonomous system
        J_res = [dg_ds,  dg_dmu];                                       % Build the Jacobian matrix
    
    elseif n_auto == 1                                  % Autonomous system 
        f_end = Fcn(0,Z_end,param);                                     % Evaluate RHS at end points of integration
        dg_domega = - 2*pi / (n_shoot * omega^2) .* reshape(f_end,dim*n_shoot,1);                   % dg/domega
        dpc_ds = f_p.';                                                 % dpc/ds (pc: Poincare phase condition, see above)
        dpc_domega = 0;                                                 % dpc/domega = 0 (phase condition is independent of the frequency)
        f_mu_plus = reshape(Fcn(0,reshape(s_p,dim,n_shoot),param_mu_plus),dim*n_shoot,1);           % Evaluate the RHS at predictor point s_p with perturbed mu
        if calc_stability                                               % Use central finite difference
            f_mu_minus = reshape(Fcn(0,reshape(s_p,dim,n_shoot),param_mu_minus),dim*n_shoot,1);     % Evaluate the RHS at predictor point s_p with perturbed mu   
            dpc_dmu = (f_mu_plus - f_mu_minus).' * (s - s_p) / (2*delta_mu);                        % dpc/dmu
        else                                                            % Use forward finite difference
            dpc_dmu = (f_mu_plus - f_p).' * (s - s_p) / delta_mu;                                   % dpc/dmu
        end
        
        J_res = [dg_ds,  dg_domega,  dg_dmu;                            % Build the Jacobian matrix
                 dpc_ds, dpc_domega, dpc_dmu];                          % Build the Jacobian matrix

    end
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