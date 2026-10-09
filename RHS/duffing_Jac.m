% This is the Jacobian of the right-hand side of the non-autonomous (forced) duffing oscillator
% The function must be able to accept a [dim x n] matrix for z, where n is the number of state vectors
% t must be a [1 x n] row vector corresponding to the state vectors z
% The Jacobian(s) must be returned as [dim x dim x n] array, where df_dz(:,:,i) is the Jacobian at z(:,i)

function df_dz = duffing_Jac(t,z,param)

    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    kappa = param{1};           % "kappa" is the first element of the "param" array
    D = param{2};               % "D" is the second element of the "param" array
    eta = param{3};             % "eta" is the third element of the "param" array
    g = param{4};               % "g" is the fourth element of the "param" array
    c = param{5};               % "c" is the fifth element of the "param" array

    n = numel(z1);              % Number of state vectors
    df_dz = NaN(2,2,n);         % Initialise
    
    % Jacobian of the RHS of the Duffing equation:
    df_dz(1,1,:) = zeros(1,n);                          % df1/dz1 = d(dz1/dt)/dz1
    df_dz(1,2,:) = ones(1,n);                           % df1/dz2 = d(dz1/dt)/dz2
    df_dz(2,1,:) = - c.*ones(1,n) - 3*kappa.*z1.^2;     % df2/dz1 = d(dz2/dt)/dz1
    df_dz(2,2,:) = - 2*D.*ones(1,n);                    % df2/dz2 = d(dz2/dt)/dz2

end