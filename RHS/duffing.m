function f = duffing(t,z,kappa,D,eta)

f(1,1) = z(2,1);
f(2,1) = -z(1,1)-kappa.*z(1,1).^3-2.*D.*z(2,1)+cos(eta.*t);

end