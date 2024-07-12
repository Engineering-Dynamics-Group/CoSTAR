% This is the right-hand side of the quasi-periodically forced harmonic oscillator

function dzdt = ho_ap_qp(t,z,param)

    % State variables and time variables
    z1 = z(1,:);
    z2 = z(2,:);
    t1 = t(1,:);
    t2 = t(2,:);

    % Parameters
    D    = param{1};
    f1   = param{2};
    f2   = param{3};
    eta = param{4};
    ratio = param{5};

    % System
    dzdt(1,:) = z2;
    dzdt(2,:) = - 2.*D.*z2 - z1 + f1.*cos(eta.*t1) + f2.*cos(ratio.*eta.*t2);

end