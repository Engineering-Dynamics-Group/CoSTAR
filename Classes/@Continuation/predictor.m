% Function predictor calculates a predictor point from the last curve point
% by using the direction_vector to the curve, which has been calculated beforehand
%
% @obj:  Continuation class object

function obj = predictor(obj)   % Compute the predictor


% Tangent or secant predictor (also used instead of parable or cubic predictor if there are less than 3 curve points)
if strcmpi(obj.pred,'tangent') || strcmpi(obj.pred,'secant') || (obj.p_local_cont_counter <= 2)

    obj.yp = obj.y0 + obj.direction.*obj.step_width.*obj.dy0;           % Calculate predictor point


% Parable predictor (also used instead of cubic predictor of there are less than 4 curve points)
elseif strcmpi(obj.pred,'parable') || (obj.p_local_cont_counter <= 3)   % Only possible if there are at least 3 curve points

    % Calculate the predictor point yp = y(sp) using the polynomial y(s) = y0 + C1*s + C2*s^2, where s is the "local" arc-length with y(s=0) = y0
    % y(-s1) = y_{k-1}  and  y(-s2) = y_{k-2}  ->  Y = [Delta_y_{k-1}, Delta_y_{k-2}]  using  Delta_y_{k-i} = y_{k-i} - y0
    % ->  Y = [C1, C2] * S  with  S = [-s1, -s2; (-s1)^2, (-s2)^2]  ->  [C1, C2] = Y/S

    % Calculate the local arc-lengths si to the last 2 solutions (with respect to y0). Why do we not use the get-function of Continuation? 
    % This is because we would have to save the y(si) AND the arc-lengths si in a Continuation property. This code here only needs to save the y(si)
    s1 = obj.direction.*norm(obj.y0 - obj.p_y0_old{3});                 % Arc-length difference between y_{k-1} = p_y0_old{3} and y_k = y0
    s2 = s1 + obj.direction.*norm(obj.p_y0_old{3} - obj.p_y0_old{2});   % Arc-length difference between y_{k-2} = p_y0_old{2} and y_k = y0
    
    % Solve the equation system for C1 and C2
    Y = [obj.p_y0_old{3}-obj.y0, obj.p_y0_old{2}-obj.y0];               % Function values for equation system
    C1 = 1/(-s1*s2^2+s2*s1^2) .* ( s2^2.*Y(:,1) - s1^2.*Y(:,2) );       % This is faster than Y/S
    C2 = 1/(-s1*s2^2+s2*s1^2) .* ( s2  .*Y(:,1) - s1  .*Y(:,2) );       % This is faster than Y/S

    % Calculate the predictor point using the parable equation yp = y(sp) = y0 + C1*sp + C2*(sp)^2
    sp = obj.direction.*obj.step_width;                                 % Arc-length at predictor point
    obj.yp = obj.y0 + C1.*sp + C2.*(sp)^2;                              % Predictor point


% Cubic predictor
elseif strcmpi(obj.pred,'cubic')                                        % Only possible if there are at least 4 curve points

    % Calculate the predictor point yp = y(sp) using the polynomial y(s) = y0 + C1*s + C2*s^2 + C3*s^3, where s is the "local" arc length with y(s=0) = y0
    % y(-si) = y_{k-i} with i \in {1,2,3}  ->  Y = [Delta_y_{k-1}, Delta_y_{k-2}, Delta_y_{k-3}]  using  Delta_y_{k-i} = y_{k-i} - y0
    % ->  Y = C * S  with  S = [-Si; (-Si).^2; (-Si).^3]  and Si = [s1, s2, s3]  ->  C = Y/S

    % Calculate the local arc-lengths si to the last 3 solutions (with respect to y0). Why do we not use the get-function of Continuation? 
    % This is because we would have to save the y(si) AND the arc-lengths si in a Continuation property. This code here only needs to save the y(si)
    s1 = obj.direction.*norm(obj.y0 - obj.p_y0_old{3});                 % Arc-length difference between y_{k-1} = p_y0_old{3} and y_k = y0
    s2 = s1 + obj.direction.*norm(obj.p_y0_old{3} - obj.p_y0_old{2});   % Arc-length difference between y_{k-2} = p_y0_old{2} and y_k = y0
    s3 = s2 + obj.direction.*norm(obj.p_y0_old{2} - obj.p_y0_old{1});   % Arc-length difference between y_{k-3} = p_y0_old{1} and y_k = y0

    % Solve the equation system
    S = [-s1, -s2, -s3; (-s1)^2, (-s2)^2, (-s3)^2; (-s1)^3, (-s2)^3, (-s3)^3];          % Matrix of arc-lengths ("x"-values)
    Y = [obj.p_y0_old{3}-obj.y0, obj.p_y0_old{2}-obj.y0, obj.p_y0_old{1}-obj.y0];       % Function values for equation system
    C = Y/S;                                                                            % Solve the equation system

    % Calculate the predictor point using the cubic equation yp = y(sp) = y0 + C1*sp + C2*(sp)^2 + C3*(sp)^3
    sp = obj.direction.*obj.step_width;                                 % Arc-length at predictor point
    obj.yp = obj.y0 + C(:,1).*sp + C(:,2).*(sp)^2 + C(:,3).*(sp)^3;     % Predictor point
    

end

end