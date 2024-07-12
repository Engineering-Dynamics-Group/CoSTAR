% Function calculates the vector of the predictor direction for the next continuation step.
% 
% @obj:     continuation class object

function obj = direction_vector(obj)
    

% Tangent predictor (also used for other predictors if there is only one curve point)
if strcmpi(obj.pred,'tangent') || (obj.p_local_cont_counter == 1)             
    [Q,~] = qr(obj.p_J0(1:end-1,:).');                                  % Do QR factorization of Jacobian without subspace constraint
    obj.dy0 = Q(:,end);                                                 % Tangent is last column of Q
      

% Secant predictor (also used for parable or cubic predictor if there are only two curve points)
elseif strcmpi(obj.pred,'secant') || (obj.p_local_cont_counter == 2)        
    v = obj.y0 - obj.p_y0_old{3};                                     % Vector from y_{k-1} = p_y0_old{3} to y_k = y0
    obj.dy0 = 1./norm(v).*v;


% Parable predictor (also used for cubic predictor if there are only three curve points)
elseif strcmpi(obj.pred,'parable') || (obj.p_local_cont_counter == 3)   
    % Calculate the coefficients of the parable y(s) = y0 + C1*s + C2*s^2, where s is the "local" arc length with y(s=0) = y0
    % -> S = [-s1, -s2; (-s1)^2, (-s2)^2];  ->  Y = [obj.p_y0_old{3}-obj.y0, obj.p_y0_old{2}-obj.y0];  ->  obj.p_C_parable = Y/S;
    s1 = obj.direction.*norm(obj.y0 - obj.p_y0_old{3});                 % Arc-length difference between y_{k-1} = p_y0_old{3} and y_k = y0
    s2 = s1 + obj.direction.*norm(obj.p_y0_old{3} - obj.p_y0_old{2});   % Arc-length difference between y_{k-2} = p_y0_old{2} and y_k = y0
    C1 = 1/(-s1*s2^2+s2*s1^2) .* ( s2^2.*(obj.p_y0_old{3}-obj.y0) - s1^2.*(obj.p_y0_old{2}-obj.y0) );     % This is faster than Y/S
    C2 = 1/(-s1*s2^2+s2*s1^2) .* ( s2  .*(obj.p_y0_old{3}-obj.y0) - s1  .*(obj.p_y0_old{2}-obj.y0) );     % This is faster than Y/S
    obj.p_C = [C1, C2];

    % The direction vector is a unit vector from y0 to yp. Therefore: Make a guess for the predictor point
    % It is only a guess since the step width still needs to be adapted (and the step control needs the direction vector)!
    % Possible solution: Take secant as direction vector when using parable predictor
    % However, this is not as accurate and lines 20 - 22 in continuation method "predictor" must be deactivated (does not work with step control)
    % obj.dy0 = (obj.y0 - obj.p_y0_old{3}) / norm(obj.y0 - obj.p_y0_old{3});   % Use secant as direction vector            
    sp_pre = obj.direction.*obj.step_width; 
    yp_pre = obj.y0 + obj.p_C(:,1).*sp_pre + obj.p_C(:,2).*(sp_pre)^2;      % Guessed predictor point
    
    % Calculate the direction vector
    obj.dy0 = (yp_pre - obj.y0) / norm(yp_pre - obj.y0);


% Cubic predictor
elseif  strcmpi(obj.pred,'cubic')                                       
    % Calculate the coefficients of the polynom y(s) = y0 + C1*s + C2*s^2 + C3*s^3, where s is the "local" arc length with y(s=0) = y0
    s1 = obj.direction.*norm(obj.y0 - obj.p_y0_old{3});                 % Arc-length difference between y_{k-1} = p_y0_old{3} and y_k = y0
    s2 = s1 + obj.direction.*norm(obj.p_y0_old{3} - obj.p_y0_old{2});   % Arc-length difference between y_{k-2} = p_y0_old{2} and y_k = y0
    s3 = s2 + obj.direction.*norm(obj.p_y0_old{2} - obj.p_y0_old{1});   % Arc-length difference between y_{k-3} = p_y0_old{1} and y_k = y0
    S = [-s1, -s2, -s3; (-s1)^2, (-s2)^2, (-s3)^2; (-s1)^3, (-s2)^3, (-s3)^3];
    Y = [obj.p_y0_old{3}-obj.y0, obj.p_y0_old{2}-obj.y0, obj.p_y0_old{1}-obj.y0];
    obj.p_C = Y/S;

    % Same problem as when using parable predictor (see lines 30 - 34 above) 
    % -> When using secant as direction vector: deactivate lines 30 - 32 in continuation method "predictor"      
    sp_pre = obj.direction.*obj.step_width; 
    yp_pre = obj.y0 + obj.p_C(:,1).*sp_pre + obj.p_C(:,2).*(sp_pre)^2 + obj.p_C(:,3).*(sp_pre)^3;      % Guessed predictor point
    
    % Calculate the direction vector
    obj.dy0 = (yp_pre - obj.y0) / norm(yp_pre - obj.y0);

end


% If direction of continuation switches between two consecutive curve points, change the direction
if (obj.p_local_cont_counter >= 2) && (obj.p_dy_old'*obj.dy0 < 0)       
    obj.dy0 = - obj.dy0;

% If direction of continuation in the initial point is negative change the sign to positiv
elseif (obj.p_local_cont_counter == 1) && (obj.dy0(end,1) < 0)          
    obj.dy0 = - obj.dy0;

end


end
