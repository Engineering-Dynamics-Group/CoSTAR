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

function   [multipliers,vectors,n_unstable,stability_flag] = PS_SHM_calc_stability_non_auto(obj,y,J,DYN,AM)
       

           stability_flag = 1;
           if strcmpi(DYN.approx_method,'shooting')  % Solution has already been computed by shooting - no need to do it again
                

                M = J(1:end-1,1:end-1) + eye(numel(y)-1);             %end-1: Exclude the row and column associated to the subspace contstraint
            

           else %Recompute the solution with shooting algorithm to get the monodromy matrix 

               %Get a starting vector in state space for the shooting method
               x0 = AM.getIC(y,DYN);   

               mu = y(end);
                %We are not solving the problem with a subspace constraint due to the following reasons:
               % 1) Any other constraint except for the natural condition (y(end)-mu0) would no guarantee meeting the value of mu0.
               % 2) Even the natural constraint does not exactly guarantee that, since it is only solved approximately. There are cases
               %    where fsolve abborts the solution process even without exactly meeting the natural constraint. This occurs when 
               %    the difference y(end)-mu0 becomes very small, which is likely, when iterating a fold bifurcation point.
               
               J = eye(numel(x0));
               try
                    [~,~,stability_flag,~,J] = fsolve(@(x)obj.SHM_single_fun(x,mu,DYN),x0,obj.fsolve_opts);      %Reshoot the solution to get J
               catch
                    stability_flag = 0;
               end
               M = J + eye(numel(x0));                   %Construct the monodromy matrix

           end

    
           [eigenvectors,eigenvalues] = eig(M);                %calcualte floquet multipliers
           [multipliers,vectors] = obj.sort_multipliers(DYN,diag(eigenvalues),eigenvectors);     %Sorting the multipliers according to different criteria
 

           cm = obj.crit_multi(DYN,multipliers);                    %Calculate the critical multipliers: This excludes automatically the Floquet multiplier = 1 for the autonomous case
           n_unstable = numel(find(cm>10*eps));                     %Number of unstable values -> 10*eps since eigenvalues on imaginary axis can have ...
                                                                    %real part which is only numerically zero but not exact

           [multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag); %Checks of NaN or Inf values 
            
     


end
 









