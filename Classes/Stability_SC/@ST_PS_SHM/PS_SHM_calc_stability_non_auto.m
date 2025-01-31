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
% @stability_flag:   Flag inidicating success of stability computation

function   [multipliers,vectors,n_unstable,stability_flag] = PS_SHM_calc_stability_non_auto(obj,y,J,DYN,AM)

dim = DYN.dim;
stability_flag = 1;

if strcmpi(DYN.approx_method,'shooting')    % Solution has already been computed by shooting - no need to do it again

    n_shoot = AM.n_shoot;

else                                        % Recompute the solution with shooting algorithm to get the monodromy matrix

    % We are not solving the problem with a subspace constraint due to the following reasons:
    % 1) Any other constraint except for the natural condition (y(end)-mu0) would no guarantee meeting the value of mu0.
    % 2) Even the natural constraint does not exactly guarantee that, since it is only solved approximately. There are cases
    %    where fsolve aborts the solution process even without exactly meeting the natural constraint. This occurs when
    %    the difference y(end)-mu0 becomes very small, which is likely, when iterating a fold bifurcation point.

    n_shoot = obj.n_shoot;
    
    % Get a starting vector in state space for the shooting method
    x0 = AM.getIC(y,DYN,n_shoot);
    mu = y(end);
    J = eye(numel(x0));                     % Preallocate

    try
        [~,~,stability_flag,~,J] = fsolve(@(x)obj.SHM_fun(x,mu,DYN),x0,obj.fsolve_opts);      % Reshoot the solution to get J
    catch
        stability_flag = 0;
    end

end

% Calculate the monodromy matrix
if n_shoot == 1
    M = J(1:dim,1:dim) + eye(dim);
else
    M = eye(dim,dim);
    for k = 1:n_shoot
        M = J((k-1)*dim+1:k*dim,(k-1)*dim+1:k*dim)*M;
    end
end

[eigenvectors,eigenvalues] = eig(full(M));      % Calculate Floquet multipliers
[multipliers,vectors] = obj.sort_multipliers(DYN,diag(eigenvalues),eigenvectors);   % Sorting the multipliers according to different criteria


cm = obj.crit_multi(DYN,multipliers);           % Calculate the critical multipliers: This excludes automatically the Floquet multiplier = 1 for the autonomous case
n_unstable = numel(find(cm>10*eps));            % Number of unstable values -> 10*eps since eigenvalues on imaginary axis can have ...
                                                % real part which is only numerically zero but not exact

[multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag);   % Checks for NaN or Inf values

end