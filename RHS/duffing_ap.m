% This is the right-hand side of the non-autonomous (forced) duffing oscillator 

function f = duffing_ap(t,z,param)

    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    kappa = param{1};           % "kappa" is the first element of the "param" array
    D = param{2};               % "D" is the second element of the "param" array
    eta = param{3};             % "eta" is the third element of the "param" array
    g = param{4};               % "g" is the fourth element of the "param" array
   
    % Duffing equation, rewritten as system of two differential equations of first order:
    f(1,:) =  z2;                                                       % First row is: d(z1)/dt = z2
    f(2,:) = - 2*D.*z2 - z1 - kappa.*z1.^3  + g.*cos(eta.*t);           % Second row is: d(z2)/dt = ...

end