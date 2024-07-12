function f = duffing_qp_test(t,z,param)

D = param{1};
kappa = param{2};
s1 = param{3};
s2 = param{4};
eta1 = param{5};
eta2 = param{6};

f(1,:) = z(2,:);
f(2,:) = -z(1,:)-kappa.*z(1,:).^3-2.*D.*z(2,:)+s1.*sin(eta1.*t)+s2.*sin(eta1.*eta2.*t);

end