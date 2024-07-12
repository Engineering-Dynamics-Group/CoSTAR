% W.F. Langford system example: Marx, Vogt - Dynamische Systeme - Theorie und Numerik p. 340
% This system shows periodic doubling bifurcations for
%omega = 3.5, rho = 0.25, beta = 0.7
%
%at the intervals epsilon \in [0.03,0.06]
% and at epsilin \in [0.0675,0.07]

function dxdt = langford(t,x,param)

    beta = param{1,1};
    omega = param{1,2};
    rho   = param{1,3};
    epsilon = param{1,4};


    dxdt = [(x(3,:)-beta).*x(1,:)-omega.*x(2,:);
             omega.*x(1,:) + (x(3,:)-beta).*x(2,:);
             0.6+x(3,:)-1./3.*x(3,:).^3-(x(1,:).^2+x(2,:).^2).*(1+rho.*x(3,:))+epsilon.*x(3,:).*x(1,:).^3];



end