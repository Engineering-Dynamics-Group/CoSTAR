% @checkReturn(y,x,iterate,exitflag) checks if the next curve point of the
% continuation approches the curve previously continued and terminates the
% continuation if this is the case
%
% Input:    DYN     DynamicalSystem class object
%           S       Solution class object
%
% Output:   obj     Continuation class object
%
function obj = check_return(obj,DYN,S)

s_mat = S.s;                                                                % Get the method solution vectors
if iscell(s_mat)                                                            % If the FGM is used, the method solution vectors are stored in cell arrays, ...
    s_mat = cell2mat(s_mat);                                                % but we need a matrix
end
if DYN.n_auto ~= 0                                                          % This if condition is need to exclude EQ from the call S.freq below
    freq_vec = S.freq(DYN.n_freq-DYN.n_auto+1:end,:);                       % Get autonomous frequencies
else
    freq_vec = [];                                                          % No frequencies for EQ
end
Y = [s_mat; freq_vec; S.mu];                                                % This is a matrix of all computed solution vectors y

N = size(Y,2);                                                              % Number of already continued curve points
y = obj.p_y1;                                                               % Current curve point

for k = 1:(N-1)
    ab = Y(:,k+1) - Y(:,k);                                                 % Vector between two consecutive curve points 
    ap = y - Y(:,k);                                                        % Vector between curve points and current curve point

    t = (ap.'*ab)/(ab.'*ab);                                                % Projection of ap onto ab
    t_clamped = max(0, min(1, t));                                          % Clipping to parameter in [0,1]

    nearest = Y(:, k) + t_clamped * ab;                                     % Next point on segment
    dist = norm(y - nearest);                                               % Euclidean distance

    if dist <= obj.p_returnTol                                              % Set error flag and messages and return if tolerance is not maintained 
        obj.p_error_flag = 0;
        obj.p_error_msg = append('ERROR: The distance of the solution for Iter = ',num2str(obj.p_local_cont_counter+1),' to the solution ', ...
                                 'curve between solutions Iter = ',num2str(k),' and Iter = ',num2str(k+1),' is <= ',num2str(obj.p_returnTol,'%.0e'),'!');
        obj.p_stopping_flag = 'CoSTAR stopped because the solution curve intersects itself.';
        return;
    end
end

end