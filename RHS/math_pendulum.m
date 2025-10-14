% This is the right-hand side of the freely oscillating mathematical pendulum

function f = math_pendulum(t,z,param)

    z1 = z(1,:);                % z1 is the first state variable
    z2 = z(2,:);                % z2 is the second state variable
                                % IMPORTANT: The state variables z_i ALWAYS have to be defined by z_i = z(i,:), e.g. z_2 = z(2,:)

    g = param{1};               % "g" is the first element of the "param" array
    l = param{2};               % "l" is the second element of the "param" array
    D = param{3};               % "D" is the third element of the "param" array
   
    % Mathematical pendulum equation, rewritten as system of two differential equations of first order:
    f(1,:) =  z2;                               % First row is: d(z1)/dt = z2
    f(2,:) = - 2*D.*z2 - g/l.*sin(z1);          % Second row is: d(z2)/dt = ...

end