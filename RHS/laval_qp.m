% This is the right-hand side of the Jeffcott(-Laval) rotor system

function f = laval_qp(t,z,param)

    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
    z3 = z(3,:);                % z3 is the third state variable
    z4 = z(4,:);                % z4 is the fourth state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)
    t1 = t(1,:);                % t1 is the first row of the time vector
  % t2 = t(2,:);                % t2 is the second row of the time vector
                                % IMPORTANT: The definition of t1 above is required for calculating quasi-periodic solutions of ...
                                %            systems exhibiting one excitation frequency!

    eta   = param{1};           % "eta" is the first element of the "param" array
    Di    = param{2};           % "Di" is the second element of the "param" array
    delta = param{3};           % "delta" is the third element of the "param" array
    e     = param{4};           % "e" is the fourth element of the "param" array
    d3    = param{5};           % "d3" is the fifth element of the "param" array
    Fg    = param{6};           % "Fg" is the sixth element of the "param" array

    % Jeffcott(-Laval) rotor system, rewritten as system of four differential equations of first order:
    f(1,:) = z3;                % First row is:  d(z1)/dt = z3
    f(2,:) = z4;                % Second row is: d(z2)/dt = z4
    f(3,:) = - 2*Di*(1+delta).*z3 - d3.*(z3.^3) - z1 + 2*Di*eta.*z2 + e*(eta^2).*sin(eta.*t1);
    f(4,:) = - 2*Di*(1+delta).*z4 - d3.*(z4.^3) - z2 - 2*Di*eta.*z1 + e*(eta^2).*cos(eta.*t1) - Fg;

    % IMPORTANT: When calculating quasi-periodic solutions, the explicit time variables belonging to different frequencies have to ...
    %            be independent! This is due to the underlaying hyper-time parametrization. 
    %            In this case, the excitation frequency is Omega_1 = eta. Thus, the time variable belonging to Omega_1 must be t1. ...
    %            t2 is not needed here since there is only one excitation frequency. t1 (and t2) ALWAYS have to be defined as above!

end