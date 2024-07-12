function f = f_coupledvdp(z,theta1,theta2,eps,alpha,beta)
   

f = eps.*[  z(1,:).^2.*z(3,:);
            z(2,:).^2.*z(4,:)];
end