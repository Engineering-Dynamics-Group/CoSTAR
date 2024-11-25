%   This is a method of the Stability subclass ST_EQ
%   It calculates the multipliers and indicates if the solution is stable or unstable        
%
%@obj:  ST_EQ object
%@y:    current solution point
%@J:    Jacobian of the current solution point
%@DYN:  Dynamical System class object
%
%@multipliers:      eigenvalues of Jacobian
%@vectors:          eigenvectors of Jacobian
%@n_unstable:       number of unstable eigenvalues
%@stability_flag:   Flag indicating success of stability computation

function   [multipliers,vectors,n_unstable,stability_flag] = EQ_calc_stability(obj,y,J,DYN)

            stability_flag = 1;     %Flag indicating success of stability computation
            tol = 10*eps;           %Tolerance to determine if value is numerically zero

            try
                [eigenvectors,eigenvalues] = eig(full(J(1:end-1,1:end-1)));     %calculate eigenvalues: The last row and column correspond to the ...
                                                                                % subspace constraint and are not relevant for the stability investigation.
                [multipliers,vectors] = obj.sort_multipliers(DYN,diag(eigenvalues),eigenvectors);     %Sorting the multipliers and corresponding vectors according to different criteria
                cm = obj.crit_multi(DYN,multipliers);                           %get the real values of the multipliers
                n_unstable = numel(find(cm>tol));                               %Number of unstable values -> tol since eigenvalues on imaginary axis can have ...
                                                                                %real part which is only numerically zero but not exact
                % if ~isempty(find(abs(real(cm))<tol,1))
                %     stability_flag = 2;         %Eigenvalue(s) on imaginary axis detected!
                % end
            catch
                multipliers = NaN(numel(y)-1,1);
                vectors = NaN(numel(y)-1);
                n_unstable = NaN;
                stability_flag = 0;
            end


           [multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag); %Checks of NaN or Inf values 
           
end
 