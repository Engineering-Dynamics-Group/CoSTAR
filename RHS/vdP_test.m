%This is the right-hand side of the autonomous (unforced) van-der-Pol
%oscillator 

function f = vdP_test(z,theta,epsilon)

%     epsilon = param(1);

    f(1,:) = -epsilon.*(z(1,:).^2-1).*z(2,:);

end