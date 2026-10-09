%% Function computing the Jacobian of Fcn via central finite differences
%
% @Fcn: RHS of the system
% @t:   Time 
% @Z:   State space vector or matrix of state space vectors [z_1, z_2, ..., z_n]
% @J:   Jacobian of Fcn with size [dim x dim x n]

function J = drhs_dz(obj,t,Z,param)
    
    % Parameters
    Fcn = obj.rhs;          % RHS of the system
    dim = size(Z,1);        % State space dimension
    n = size(Z,2);          % Number of state vectors
    
    h = eps^(1/3);
    Z_vec = reshape(Z,n*dim,1);
    
    % Set up the perturbed z_i vectors: Since each z_i needs to be perturbed dim-times, we get dim pertubed vectors for z_i,
    % where the j-th vector is the perturbed z_i in the j-th component. This is done for all n state vectors
    % The H matrix stores the individual step widths h_(i,j) = h*(1+abs(z_(i,j)) for numerical differentiation in a [dim x dim*n] matrix
    Z_dim = reshape(repmat(Z,dim,1),dim,n*dim);         % This is a [dim x dim*n] matrix where each z_i is repeated dim times
    H = h.*(repmat(eye(dim),1,n) + sparse(repmat(1:1:dim,1,n),1:1:n*dim,abs(Z_vec),dim,n*dim));
    Z_dim_plus  = Z_dim + H;                            % This is a [dim x dim*n] matrix containing all "+ H" perturbed z_i
    Z_dim_minus = Z_dim - H;                            % This is a [dim x dim*n] matrix containing all "- H" perturbed z_i
    t_dim = reshape(repmat(t,dim,1),1,n*dim);           % Since each z_i occurs dim times, we also need each t_i value dim times
    
    % Compute the Jacobian vectorized for all n state vectors. First, J is a [dim x dim*n] matrix, where the n Jacobian(s) are placed next to each other
    if strcmpi(obj.sol_type,'equilibrium')              % The RHS of equilibria lacks the dependence on t, which is why the call of Fcn is slightly different
        J = (Fcn(Z_dim_plus,param) - Fcn(Z_dim_minus,param)) ./ (2.*nonzeros(H).');
    else
        J = (Fcn(t_dim,Z_dim_plus,param) - Fcn(t_dim,Z_dim_minus,param)) ./ (2.*nonzeros(H).');
    end
    J = reshape(J,dim,dim,n);                           % Reshape J to the desired output structure
    
end