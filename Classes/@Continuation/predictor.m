% Function predictor calculates a predictor point from the last curve point
% by using the direction_vector to the curve, which has been calculated beforehand
%
%@obj:  Continuation class object

function obj = predictor(obj)   % Compute the predictor


% Tangent or secant predictor (also used instead of parable or cubic predictor if there are less than 3 curve points)
if strcmpi(obj.pred,'tangent') || strcmpi(obj.pred,'secant') || (obj.p_local_cont_counter <= 2)

    obj.yp = obj.y0 + obj.direction.*obj.step_width.*obj.dy0;           % Calulate predictor point


% Parable predictor (also used instead of cubic predictor of there are less than 4 curve points)
elseif strcmpi(obj.pred,'parable') || (obj.p_local_cont_counter <= 3)   % Only possible if there are at least 3 curve points

    sp = obj.direction.*obj.step_width;
    obj.yp = obj.y0 + obj.p_C(:,1).*sp + obj.p_C(:,2).*(sp)^2;          % Calculate acutal predictor point after step width has been adapted
    obj.dy0 = (obj.yp - obj.y0) / norm(obj.yp - obj.y0);                % Calculate the correct direction vector after the predictor point has been determined
    % Evtl. insgesamt verbessern -> dy0 muss hier nicht berechnet werden

% Cubic predictor
elseif strcmpi(obj.pred,'cubic')                                        % Only possible if there are at least 4 curve points

    sp = obj.direction.*obj.step_width;
    obj.yp = obj.y0 + obj.p_C(:,1).*sp + obj.p_C(:,2).*(sp)^2 + obj.p_C(:,3).*(sp)^3;  % Calculate acutal predictor point after step width has been adapted
    obj.dy0 = (obj.yp - obj.y0) / norm(obj.yp - obj.y0);                % Calculate the correct direction vector after the predictor point has been determined
    % Evtl. insgesamt verbessern -> dy0 muss hier nicht berechnet werden
    
end


end