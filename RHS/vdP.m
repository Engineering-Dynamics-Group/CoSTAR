function f = vdP(t,z,epsilon,eta)

f(1,1) = z(2,1);
f(2,1) = -z(1,1)-epsilon.*(z(1,1).^2-1).*z(2,1)+cos(eta.*t);

end

