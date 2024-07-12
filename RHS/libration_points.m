% This is the function to determine the libration points L1, L2 and L3

function f = libration_points(z,param)

    z1 = z(1,:);                    % z1 is the only state variable
                                    % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    mu = param{1};                  % "mu" is the first element of the "param" array
    point = param{2};

    % This is the equation relevant for L1, L2 and L3:
    % f = - (1-mu).*(z1+mu)./(abs(z1+mu).^3) - mu.*(z1-1+mu)./(abs(z1-1+mu).^3) + z1;
    % The desired libration point to compute is "chosen" by the initial value

    % In order to obtain a unique solution, the abs() function is replaced by s_d*() and s_r*(), where s_d and s_r equals +/- 1
    % The s_d and s_r are chosen according to the desired libration point
    % This should make sure that a different libration point can not be computed when the initial value is too far away from ...
    % the desired libration point
    switch point
        case 'L1'
            s_d = 1;    s_r = -1;
        case 'L2'
            s_d = 1;    s_r = 1;
        case 'L3'
            s_d = -1;    s_r = -1;
    end

    f = - (1-mu).*(z1+mu)./(s_d.*(z1+mu).^3) - mu.*(z1-1+mu)./(s_r.*(z1-1+mu).^3) + z1;     % This is the modified equation
    
end