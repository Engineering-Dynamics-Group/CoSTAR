function Z0 = initializeShootingNodes(obj,Z_C0,FW,T)
%INITIALIZESHOOTINGNODES
% Generates an initial guess for all multiple-shooting nodes from the
% initial characteristic states Z_C0 on the boundary.
% Returns the complete [dim, n_shoot, n_char]-matrix that is needed for
% the (multiple-)shooting method.
%
% Z_C0(:,k)   = initial state of characteristic k
% Z0(:,i,k)   = initial state of sub-interval i of characteristic k
%
% Dimensions:
%   Z_C0 : [dim x n_char]
%   Z0   : [dim x n_shoot x n_char]

%% Parameters

dim = obj.n;
n_char = obj.n_char;
n_shoot = obj.n_shoot;

%% Check dimensions

if ~isequal(size(Z_C0),[dim,n_char])
    error(['Z_C0 must have size [dim,nChar] = [', ...
           num2str(dim),',',num2str(n_char),'].']);
end

%% Initialize complete initial-value matrix

Z0 = zeros(dim,n_shoot,n_char);
Z0(:,1,:) = reshape(Z_C0,[dim,1,n_char]);                                   % Embed the given initial values on the boundary into Z0

%% No internal nodes required for single shooting

if n_shoot == 1
    return
end

%% Divide the characteristic into sub-intervals

tSubIncr = T/n_shoot;                                                       % Length of a sub-interval
tSub1 = 0;                                                                  % Start value for first sub-interval
tSub2 = tSubIncr;                                                           % End value for first sub-interval

Z_current = Z_C0;                                                           % Given IVs on the boundary are used as first IVs for the forward integration 

%% Generate the internal shooting-node guesses via forward integration

for i = 2:n_shoot

    [~,char_trajectory] = obj.solver_function( ...                          % Solves the ODE with RHS FW and initial values Z_current(:)
        FW,[tSub1,tSub2],Z_current(:),obj.odeOpts);                         % on the interval [tSub1,tSub2] 

    Z_current = reshape(char_trajectory(end,:).',[dim,n_char]);             % Set the calculated point as the new IV for the forward integration

    Z0(:,i,:) = reshape(Z_current,[dim,1,n_char]);                          % Embed the guess into Z0

    tSub1 = tSub2;                                                          % Set start value for new sub-interval
    tSub2 = tSub2 + tSubIncr;                                               % Set end value for new sub-interval

end

end