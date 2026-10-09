% This is the Jacobian of the right-hand side of the autonomous (unforced) van der Pol oscillator
% The function must be able to accept a [dim x n] matrix for z, where n is the number of state vectors
% t must be a [1 x n] row vector corresponding to the state vectors z
% The Jacobian(s) must be returned as [dim x dim x n] array, where df_dz(:,:,i) is the Jacobian at z(:,i)

function df_dz = vdP_auto_Jac(t,z,param)

    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    epsilon = param{1};         % epsilon is the first (and only) element of the "param" array

    n = numel(z1);              % Number of state vectors
    df_dz = NaN(2,2,n);         % Initialise
    
    % Jacobian of the RHS of the van der Pol equation:
    df_dz(1,1,:) = zeros(1,n);                  % df1/dz1 = d(dz1/dt)/dz1
    df_dz(1,2,:) = ones(1,n);                   % df1/dz2 = d(dz1/dt)/dz2
    df_dz(2,1,:) = - 2*epsilon.*z1.*z2 - 1;     % df2/dz1 = d(dz2/dt)/dz1
    df_dz(2,2,:) = - epsilon.*(z1.^2 - 1);      % df2/dz2 = d(dz2/dt)/dz2

end