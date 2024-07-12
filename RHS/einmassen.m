function f = einmassen(t,z,D,eta)

f(1,1) = z(2,1);
f(2,1) = -z(1,1)-2.*D.*z(2,1)+cos(eta.*t);

end
