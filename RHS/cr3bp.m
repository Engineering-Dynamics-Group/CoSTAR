% This is the function of the right-hand side of the circular restricted three-body problem (CR3BP)
%
% The location of the third body of negligible mass is charactarised by the coordinates (x,y,z) with respect to a rotating ...
% reference frame sitting at the barycentre of the system. The reference frame rotates at constant speed omega and its basis ...
% vector e_x always lies on the connecting line of the two primary bodies P1 and P2 (so called primaries), pointing at body ...
% P2 with mass m_2 <= m_1.
% The equations of motion are completely nondimensionalised with respect to the distance between the primaries and, considering ...
% time t, with respect to 1/omega (i.e. one orbital period of the primaries equals t_T = T * omega = 2*pi). Furthermore, the ...
% mass ratio mu = m_2/(m_1+m_2) is introduced.

function f = cr3bp(t,z_vec,param)

    % State variables
    x = z_vec(1,:);         % x coordinate
    y = z_vec(2,:);         % y coordinate
    z = z_vec(3,:);         % z coordinate
    v_x = z_vec(4,:);       % velocity dx/dt
    v_y = z_vec(5,:);       % velocity dy/dt
    v_z = z_vec(6,:);       % velocity dz/dt


    % Parameter
    mu = param{1};      % Mass ratio

    d = vecnorm([x+mu; y; z]);      % Distance from primary P1 to third body
    r = vecnorm([x-1+mu; y; z]);    % Disctance from primary P2 to third body
    % Note: vecnorm automatically computes the Euclidean norm of each column of a matrix (the position "vectors" can be matrices, see definition above)


    % First-order system
    f(1,:) = v_x;
    f(2,:) = v_y;
    f(3,:) = v_z;
    f(4,:) = - (1-mu) ./ (d.^3) .* (x + mu)  -  mu ./ (r.^3) .* (x - 1 + mu) + x;
    f(5,:) = - (1-mu) ./ (d.^3) .* y         -  mu ./ (r.^3) .* y            + y;
    f(6,:) = - (1-mu) ./ (d.^3) .* z         -  mu ./ (r.^3) .* z;

end