%This is the right-hand side of the autonomous (unforced) van-der-Pol
%oscillator 

function f = vdP_forced(t,z,param)

    epsilon = param{1};
    s = param{2};
    eta = param{3};

    f(1,:) = z(2,:);
    f(2,:) = -z(1,:)-epsilon.*(z(1,:).^2-1).*z(2,:)+s.*cos(eta.*t);

end