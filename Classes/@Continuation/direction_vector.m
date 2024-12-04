% This Continuation method calculates a direction vector, which is used for tangent or ... 
% secant predictor (if they were chosen) and as geometrical information for step control.
% 
% @obj:     continuation class object

function obj = direction_vector(obj)
    

% Case: Tangent predictor or there is only one curve point
% if strcmpi(obj.pred,'tangent') || ((obj.p_local_cont_counter == 1)&&(strcmpi(obj.pred,'parable')||strcmpi(obj.pred,'cubic'))) || obj.p_use_qr        
if strcmpi(obj.pred,'tangent') || obj.p_use_qr        
    [Q,~] = qr(obj.p_J0(1:end-1,:).');          % Do QR factorization of Jacobian without subspace constraint
    obj.dy0 = Q(:,end);                         % Tangent is last column of Q
    obj.p_use_qr = false;                       % Reset property since secant predictor can now be used

% Case: All other predictors (secant, parable, cubic) and there is only 1 curve point
elseif obj.p_local_cont_counter == 1 
    v = obj.p_initial_slope;                    % Vector calculated to determine initial slope by secant
    obj.dy0 = 1./norm(v).*v;

% Case: All other predictors (secant, parable, cubic) and there are at least 2 curve points 
else
    v = obj.y0 - obj.p_y0_old{3};               % Vector from y_{k-1} = p_y0_old{3} to y_k = y0
    obj.dy0 = 1./norm(v).*v;

end


% If the direction of continuation switches between two consecutive curve points: change the direction
if (obj.p_local_cont_counter >= 2) && (obj.p_dy_old'*obj.dy0 < 0)       
    obj.dy0 = - obj.dy0;


% If the direction of continuation at the initial point is negative: change the sign to positive
elseif (obj.p_local_cont_counter == 1) && (obj.dy0(end,1) < 0)          
    obj.dy0 = - obj.dy0;


end


end