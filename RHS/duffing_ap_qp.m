% This is the right-hand side of the quasi-periodically forced (non-autonomous) duffing oscillator 

function f = duffing_ap_qp(t,z,param)
    
    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)
    t1 = t(1,:);                % t1 is the first row of the time vector
    t2 = t(2,:);                % t2 is the second row of the time vector
                                % IMPORTANT: The definition of t1 and t2 above is required for calculating quasi-periodic solutions ...
                                %            of systems exhibiting two excitation frequencies!

    D     = param{1};           % "kappa" is the first element of the "param" array
    kappa = param{2};           % "D" is the second element of the "param" array
    f1    = param{3};           % "f1" is the third element of the "param" array
    f2    = param{4};           % "f2" is the forth element of the "param" array
    eta   = param{5};           % "eta" is the fifth element of the "param" array
    ratio = param{6};           % "ratio" is the sixth element of the "param" array

    % Duffing equation, rewritten as system of two differential equations of first order:
    f(1,:) =  z2;                                                                              % First row is: d(z1)/dt = z2
    f(2,:) = - 2*D.*z2 - z1 - kappa.*z1.^3 + f1.*sin(eta.*t1) + f2.*sin(eta.*ratio.*t2);       % Second row is: d(z2)/dt = ...

    % IMPORTANT: When calculating quasi-periodic solutions, the explicit time variables belonging to different frequencies have ...
    %            to be independent! This is due to the underlaying hyper-time parametrization. 
    %            In this case, the excitation frequencies are Omega_1 = eta and Omega_2 = ratio*eta. Thus, the time variable ...
    %            belonging to Omega_1 must be t1 and the time variable belonging to Omega_2 must be t2. t1 and t2 ALWAYS have to ...
    %            be defined as above!

end
