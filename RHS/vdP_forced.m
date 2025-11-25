function f = vdP_forced(t,z,param)

epsilon = param{1};
s1 = param{2};
eta1 = param{3};

f(1,:) = z(2,:);
f(2,:) = -z(1,:)-epsilon.*(z(1,:).^2-1).*z(2,:)+s1.*cos(eta1.*t(1,:));

end