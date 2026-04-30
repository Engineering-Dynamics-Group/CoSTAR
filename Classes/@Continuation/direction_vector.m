% This Continuation method calculates a direction vector, which is used for tangent or ... 
% secant predictor (if they were chosen) and as geometrical information for step control.
% 
% @obj:     continuation class object

function obj = direction_vector(obj,DYN)
    

% Case: Tangent predictor or there is only one curve point
% if strcmpi(obj.pred,'tangent') || ((obj.p_local_cont_counter == 1)&&(strcmpi(obj.pred,'parable')||strcmpi(obj.pred,'cubic'))) || obj.p_use_qr        
if strcmpi(obj.pred,'tangent') || obj.p_use_qr
    if isfield(DYN.system,'first_integral') 
        Jg = full(obj.p_J0(1:end-1,:));         % Jg needs to have a defect = 1 and its null space can be computed by means of a SVD
        [~,S,V] = svd(Jg);                      % The column of V corresponding to the SV = 0 is a tangent vector
        tol = 1e-8;   adapt_tol = true;         % However, the smallest SV is often not exactly zero, which is why we need the tolerance tol   
        while adapt_tol              
            idx_null = find(diag(S) < tol);             % Find the indices of SV smaller than tol
            defect = numel(idx_null);                   % Get the number of SV smaller than tol
            if      defect == 1;  adapt_tol = false;    % A single index was found - now we can get the tangent from V
            elseif  defect == 0;  tol = 10*tol;         % Increase the tolerance if there is no SV
            elseif  defect > 1;   tol = tol/10;         % Decrease the tolerance if there are too many SV 
            end
            if tol == 1e-5;       error('The tolerance for computing the tangent has reached 1e-5! Please take a look at the Jacobian!');
            elseif tol == 1e-17;  error('The tolerance for computing the tangent has reached 1e-17! Please take a look at the Jacobian!');
            end
        end
        obj.dy0 = V(:,idx_null);                % Get the tangent vector from V. norm(dy0) is already 1, since V is unitary
    else
        [Q,~] = qr(obj.p_J0(1:end-1,:).');      % Do QR factorization of Jacobian without subspace constraint
        obj.dy0 = Q(:,end);                     % Tangent is last column of Q
    end
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

% J = full(obj.p_J0);
% Use G and PC                        
% J1 = J(1:end-2,:).'; [Q1,R1] = qr(J1); t1 = Q1(:,end);
% Use G, PC and IC
% J2 = J(1:end-1,:).'; [Q2,R2] = qr(J2); t2 = Q2(:,end);
% Use G and IC
% J3 = [J(1:end-3,:); J(end-1,:)].'; [Q3,R3] = qr(J3); t3 = Q3(:,end);
% % Use G, IC and PC
% J4 = [J(1:end-3,:); J(end-1,:); J(end-2,:)].'; [Q4,R4] = qr(J4); t4 = Q4(:,end);
% J3 = full(obj.p_J0(1:end-2,:)); [Q3,R3] = qr(J3.'); obj.t = Q3(:,end);

J5 = full(obj.p_J0(1:end-1,:));   tol = 1e-8;   increase_tol = true;
while increase_tol
    defect = size(J5,2) - rank(J5,tol);
    if defect == 0;     tol = 10*tol;  
    elseif defect == 1; increase_tol = false;
    else;               error('Defect of J is not 0 or 1!');               
    end
    if tol == 1e-5;  error('Tolerance to compute the nullspace of J reached 1e-5!');  end
end
obj.t = null(J5,tol);
if obj.dy0.'*obj.t < 0;  obj.t = -obj.t;  end
% if (obj.p_local_cont_counter == 1)
%     obj.t = NaN(size(obj.dy0));
% else
%     v = obj.y0 - obj.p_y0_old{3};               % Vector from y_{k-1} = p_y0_old{3} to y_k = y0
%     obj.t = 1./norm(v).*v;
% end

end