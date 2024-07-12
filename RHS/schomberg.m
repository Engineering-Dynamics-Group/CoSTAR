function dzdt = schomberg(t,z,param)

    c = param{1,1};
    b = param{1,2};
    a = param{1,3};
    omega = param{1,4};


    dzdt = a./c.*abs(sin(omega.*t))-b./c.*z;



end