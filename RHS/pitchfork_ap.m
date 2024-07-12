% This is the function of the pitchfork equation

function f = pitchfork_ap(z,param)

    z1 = z(1,:);                        % z1 is the only state variable
                                        % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    mu = param{1};                      % "mu" is the first element of the "param" array
    gamma = param{2};                   % "gamma" is the second element of the "param" array

    f = mu.*z1 - z1.^3 + gamma;         % This is the equation

end