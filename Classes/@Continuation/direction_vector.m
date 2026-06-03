% This Continuation method calculates a direction vector, which is used for tangent or ... 
% secant predictor (if they were chosen) and as geometrical information for step control.
% 
% @obj:     continuation class object

function obj = direction_vector(obj,DYN)


% Case: Tangent predictor or initial_slope failed
if strcmpi(obj.pred,'tangent') || obj.p_use_qr
    if isfield(DYN.system,'first_integral') 
        Jg = full(obj.p_J0(1:end-1,:));         % Jg needs to have a defect = 1. Its null space can be computed by means of a SVD
        [~,S,V] = svd(Jg);                      % The column of V corresponding to the SV = 0 is a tangent vector
        tol = max(size(Jg)) * eps(norm(Jg));    % The smallest SV is often not exactly zero, which is why we need the tolerance tol
        adapt_tol = true;  counter = 0;         % tol above is the default tolerance used by the null or rank function
        while adapt_tol
            counter = counter + 1;                      % counter is used to stop the while loop in the unlikely case that tol is in- and decreased alternately
            idx_null = find(diag(S) < tol);             % Find the indices of SV smaller than tol
            defect = numel(idx_null);                   % Get the number of SV smaller than tol
            if      defect == 1;  adapt_tol = false;    % A single index was found - now we can get the tangent from V
            elseif  defect == 0;  tol = 10*tol;         % Increase the tolerance if there is no SV = 0
            elseif  defect > 1;   tol = tol/10;         % Decrease the tolerance if there are too many SV = 0
            end
            if counter > 10;        error('Computation of tangent via SVD of Jacobian Jg not possible, too many iterations for adapting the tolerance. Please check the Jacobian.')
            elseif tol >= 1e-5;     error('Computation of tangent via SVD of Jacobian Jg not possible, since the null space seems to be empty (tolerance for singular value is >= 1e-6)!');
            elseif tol <= eps/10;   error('Computation of tangent via SVD of Jacobian Jg not possible, since dimension of the null space seems to be > 1 (tolerance for singular value is <= eps)!');
            end
        end
        obj.dy0 = V(:,idx_null);                % Get the tangent vector from V. idx_null should be equal to size(V,2). norm(dy0) is already 1, since V is unitary
    else
        [Q,~] = qr(obj.p_J0(1:end-1,:).');      % Do QR factorization of Jacobian without subspace constraint
        obj.dy0 = Q(:,end);                     % Tangent is last column of Q
    end
    obj.p_use_qr = false;                  % Reset property since secant predictor can now be used

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