% This is the function of the parable equation

function f = parable(z,param)

    z1 = z(1,:);                    % z1 is the only state variable
                                    % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    mu = param{1};                  % "mu" is the first element of the "param" array
    a = param{2};                   % "a" is the second element of the "param" array
    b = param{3};                   % "b" is the third element of the "param" array

    f = mu - b + a.*z1.^2;          % This is the equation

end