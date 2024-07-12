% This is the right-hand side of the autonomous (unforced) van der Pol oscillator 

function f = vdP_auto_ap(t,z,param)

    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    epsilon = param{1};         % epsilon is the first (and only) element of the "param" array

    % van der Pol oscillator, rewritten as system of two differential equations of first order:
    f(1,:) = z2;
    f(2,:) = - epsilon .* (z1.^2 - 1) .* z2 - z1;

end