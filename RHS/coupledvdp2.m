function dz_dt = coupledvdp2(t,z,param)
   
lambda1 = param{1,1};
lambda2 = param{1,2};
omega1 = param{1,3};
omega2 = param{1,4};
BR = param{1,5};
BD = param{1,6};

dz_dt = [ z(3,:);
          z(4,:);
         (lambda1-z(1,:).^2).*z(3,:)-omega1.^2.*z(1,:)+BR.*(z(2,:)-z(1,:))+BD.*(z(4,:)-z(3,:));
         (lambda2-z(2,:).^2).*z(4,:)-omega2.^2.*z(2,:)+BR.*(z(1,:)-z(2,:))+BD.*(z(3,:)-z(4,:))];
end