% This is the right-hand side of the coupled van der Pol oscillator 

function f = coupledvdp(t,z,param)
   
    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
    z3 = z(3,:);                % z3 is the third state variable
    z4 = z(4,:);                % z4 is the fourth state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)
  % t1 = t(1,:);                % t1 is the first row of the time vector
  % t2 = t(2,:);                % t2 is the second row of the time vector
                                % IMPORTANT: The definition of t1 (and t2) above is required for calculating quasi-periodic solutions ...
                                %            of non-autonomous (and/or mixed) systems (i.e. there are one or two excitation frequencies)!

    epsilon = param{1,1};       % "epsilon" is the first element of the "param" array
    alpha   = param{1,2};       % "alpha" is the second element of the "param" array
    beta    = param{1,3};       % "beta" is the third element of the "param" array

    % coupled van der Pol oscillator, rewritten as system of four differential equations of first order:
    f(1,:) = z3;                                                                    % First row is:  d(z1)/dt = z3
    f(2,:) = z4;                                                                    % Second row is: d(z2)/dt = z4
    f(3,:) = - epsilon.*(z1.^2 - 1).*z3 -           z1 + alpha.*(z2 - z1);          % Third row is:  d(z3)/dt = ...
    f(4,:) = - epsilon.*(z2.^2 - 1).*z4 - (1+beta).*z2 + alpha.*(z1 - z2);          % Fourth row is: d(z4)/dt = ...

    % IMPORTANT: When calculating quasi-periodic solutions, the explicit time variables belonging to different frequencies have to ...
    %            be independent! This is due to the underlaying hyper-time parametrization. 
    %            In this case, there are no excitation frequencies, so we do not have to consider this.

end