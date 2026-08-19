% @checkReturn(y,x,iterate,exitflag) checks if the next curve point of the
% continuation approches the curve previously continued and terminates the
% continuation if this is the case
%
% Input:    y           New curve point
%           x           List of all previous curve points
%
% Output:   []
%
function [] = checkReturn(obj,x)

if(~obj.p_contDo); return; end                                              % If iterate is already set to false return and keep exitflag

tolReturn = 1e-3;                                                           % Tolerance for detecting approach to curve
N = size(x,2);                                                              % Number of already continued curve points

y = obj.p_y1;                                                               % Current curve point

for k = 1:(N - 1)
    ab = x(:,k+1) - x(:,k);                                                 % Vector between two consecutive curve points 
    ap = y - x(:,k);                                                        % Vector between curve points and current curve point

    t = (ap.'*ab)/(ab.'*ab);                                                % Projection of ap onto ab
    t_clamped = max(0, min(1, t));                                          % Clipping to parameter in [0,1]

    nearest = x(:, k) + t_clamped * ab;                                     % Next point on segment
    dist = norm(y - nearest);                                               % Euclidean distance

    if(dist<= tolReturn)                                                    % Set iterator and stopping flag and return if tolerance is not maintained 
        % obj.p_contDo = false;    exitflag = 2;
	    obj.p_contDo = false;    obj.p_stopping_flag = 'CoSTAR stopped, because the curve closed a loop.';
        return;
    end
end

end