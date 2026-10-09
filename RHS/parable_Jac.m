% This is the Jacobian of the right-hand side of the parable equation
% The function must be able to accept a [dim x n] matrix for z, where n is the number of state vectors
% t must be a [1 x n] row vector corresponding to the state vectors z
% The Jacobian(s) must be returned as [dim x dim x n] array, where df_dz(:,:,i) is the Jacobian at z(:,i)

function df_dz = parable_Jac(z,param)

    z1 = z(1,:);                    % z1 is the only state variable
                                    % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    mu = param{1};                  % "mu" is the first element of the "param" array
    a = param{2};                   % "a" is the second element of the "param" array
    b = param{3};                   % "b" is the third element of the "param" array

    df_dz = 2*a.*z1;                % Jacobian of the parable equation

end