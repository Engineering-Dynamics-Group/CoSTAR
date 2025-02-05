% This is a method of the Stability subclass ST_PS_SHM
% It calculates the Floquet multipliers and indicates if the solution is stable or unstable
%
% @obj:  ST_PS_SHM object
% @y:    current solution point
% @J:    Jacobian of the current solution point
% @DYN:  Dynamical System class object
% @AM:   Approximation Method subclass object
%
% @multipliers:      Floquet multipliers
% @vectors:          Eigenvectors corresponding to Floquet multipliers
% @n_unstable:       Number of unstable Floquet multipliers
% @stability_flag:   Flag indicating success of stability computation

function   [multipliers,vectors,n_unstable,stability_flag] = PS_SHM_calc_stability(obj,y,J,DYN,AM)

% Parameter
dim = DYN.dim;                              % Dimension of the system
n_auto = DYN.n_auto;                        % Number of autonomous frequencies
stability_flag = 1;                         % Initialize stability flag


%% Get the Jacobian

if strcmpi(DYN.approx_method,'shooting')    % Solution has already been computed by shooting - no need to do it again
    
    n_shoot = AM.n_shoot;

else    % Recompute the solution with shooting algorithm to get the monodromy matrix

    % We are not solving the problem with a subspace constraint due to the following reasons:
    % 1) Any other constraint except for the natural condition (y(end)-mu0) would no guarantee meeting the value of mu0.
    % 2) Even the natural constraint does not exactly guarantee that, since it is only solved approximately. There are cases
    %    where fsolve aborts the solution process even without exactly meeting the natural constraint. This occurs when
    %    the difference y(end)-mu0 becomes very small, which is likely when iterating a fold bifurcation point.

    n_shoot = obj.n_shoot;                      % Number of shooting points defined by stability class
    mu = y(end);                                % Continuation parameter
    
    x0 = AM.getIC(y,DYN,n_shoot);               % Get a starting vector in state space for the shooting method 
    if n_auto == 1
        x0 = [x0; y(end-1)];                    % Add the autonomous frequency if system is autonomous
    end
    
    J = eye(dim*n_shoot);                       % Preallocate Jacobian
    try
        [~,~,stability_flag,~,J] = fsolve(@(x) PS_SHM_ST_residuum(obj,x,x0,mu,DYN), x0, obj.fsolve_opts);   % Reshoot the solution to get J
    catch
        stability_flag = 0;                     % If shooting did fail for some reasons, flag is set to 0
    end
    
end


%% Calculate the monodromy matrix

if n_shoot == 1
    M = J(1:dim,1:dim) + eye(dim);
else
    M = eye(dim,dim);
    for k = 1:n_shoot
        M = J((k-1)*dim+1:k*dim,(k-1)*dim+1:k*dim)*M;
    end
end


%% Calculate Floquet multipliers, corresponding eigenvectors and number of unstable multipliers

[eigenvectors,eigenvalues] = eig(full(M));      % Calculate eigenvectors and eigenvalues of the monodromy matrix
[multipliers,vectors] = obj.sort_multipliers(DYN,diag(eigenvalues),eigenvectors);   % Sorting the multipliers according to different criteria

cm = obj.crit_multi(DYN,multipliers);           % Calculate the critical multipliers: This excludes automatically the multiplier = 1 in the autonomous case
n_unstable = numel(find(cm>10*eps));            % Number of unstable values -> 10*eps since eigenvalues on imaginary axis can have ...
                                                % real part which is only numerically zero but not exact

[multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag);   % Checks for NaN or Inf values


end