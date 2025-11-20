% This is the right-hand side of the periodically forced harmonic oscillator

function dzdt = harmonic(t,z,param)

    % State variables
    z1 = z(1,:);
    z2 = z(2,:);

    % Parameters
    D   = param{1};
    f   = param{2};
    eta = param{3};
    % D = 0.2;   f = 1;   eta = 0.1;

    % System
    dzdt(1,:) = z2;
    dzdt(2,:) = - 2.*D.*z2 - z1 + f.*cos(eta.*t);

end