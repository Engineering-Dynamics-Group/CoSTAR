function f = f_duffing_QPS(z,theta1,theta2,gamma,g)

f = gamma.*z(1,:).^3 - g.*(sin(theta1)+ sin(theta2));

end