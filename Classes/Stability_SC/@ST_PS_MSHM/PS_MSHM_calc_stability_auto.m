%   This is a method of the Stability subclass ST_PS_SHM
%   It calculates the Floquet multipliers and indicates if the solution is stable or unstable
%
%@obj:  ST_PS_SHM object
%@y:    current solution point
%@J:    Jacobian of the current solution point
%@DYN:  Dynamical System class object
%@AM:   Approximation Method subclass object
%
%@multipliers:      Floquet multipliers
%@vectors:          eigenvectors corresponding to Floquet multipliers
%@n_unstable:       number of unstable Floquet multipliers
%@stability_flag:   Flag inidicating success of stability computation

function   [multipliers,vectors,n_unstable,stability_flag] = PS_MSHM_calc_stability_auto(obj,y,J,DYN,AM)

stability_flag = 1;
if strcmpi(DYN.approx_method,'mshm')  % Solution has already been computed by shooting - no need to do it again

    dim = DYN.dim;
    n_shoot = AM.n_shoot;

    M0 = eye(dim,dim);
    for k=1:n_shoot
        M0 = J((k-1)*dim+1:k*dim,(k-1)*dim+1:k*dim)*M0;
    end
    M = M0;

else %Recompute the solution with shooting algorithm to get the monodromy matrix


    %We are not solving the problem with a subspace constraint due to the following reasons:
    % 1) Any other constraint except for the natural condition (y(end)-mu0) would no guarantee meeting the value of mu0.
    % 2) Even the natural constraint does not exactly guarantee that, since it is only solved approximately. There are cases
    %    where fsolve abborts the solution process even without exactly meeting the natural constraint. This occurs when
    %    the difference y(end)-mu0 becomes very small, which is likely, when iterating a fold bifurcation point.

    dim = DYN.dim;
    n_shoot = obj.n_shoot;
    M0 = eye(dim,dim);

    %Get a starting vector in state space for the shooting method
    x0 = [AM.getIC(y,DYN,n_shoot);y(end-1)];            %Add the autonomous frequency to the state space initial condition.
    mu = y(end);
    J = eye(numel(x0)-1);                       %Prealloquate
    try
        [~,~,stability_flag,~,J] = fsolve(@(x) MSHM_auto_fun(obj,x,x0,mu,DYN),x0);  %solve corrector-function
    catch
        stability_flag = 0;   %If shooting did fail for some reasons, flag is set to 0
    end

    for k=1:n_shoot
        M0 = J((k-1)*dim+1:k*dim,(k-1)*dim+1:k*dim)*M0;
    end
    M = M0;
end


[eigenvectors,eigenvalues] = eig(M);            %calcualte floquet multipliers
[multipliers,vectors] = obj.sort_multipliers(DYN,diag(eigenvalues),eigenvectors);     %Sorting the multipliers according to different criteria

cm = obj.crit_multi(DYN,multipliers);           %Calculate the critical multipliers: This excludes automatically the Floquet multiplier = 1 for the autonomous case
n_unstable = numel(find(cm>10*eps));            %Number of unstable values -> 10*eps since eigenvalues on imaginary axis can have ...
%real part which is only numerically zero but not exact

[multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag); %Checks of NaN or Inf values

end









