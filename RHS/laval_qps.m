function dxdt = laval_qps(t,x,param)

eta = param{1,1};
Di = param{1,2};
Delta = param{1,3};
e = param{1,4};
d3 = param{1,5};
Fg = param{1,6};


q1 = x(1,:);
q2 = x(2,:);
v1 = x(3,:);
v2 = x(4,:);


dxdt = [
    v1;
    v2;
    -q1 + 2*Di*eta.*q2 - 2*Di*(1+Delta).*v1 - d3.*v1.^3 + e*eta^2*sin(eta.*t(1,:));
    -q2 - 2*Di*eta.*q1 - 2*Di*(1+Delta).*v2 - d3.*v2.^3 + e*eta^2*cos(eta.*t(1,:)) - Fg];

end







