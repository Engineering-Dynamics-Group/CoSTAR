%When called, this method gives back a vector of the values of a test function (or multipliers) that looses its stability
%If no test function is available or if the test functions do not have a zero crossing, then the multipliers are investigated
%
%@obj:          Stability subclass object
%@DYN:          DynamicalSystem class object   
%
%@stab_fcn:     array of the stability function values


function stab_fcn = get_stability_fcn(obj,DYN)

    curve_container = obj.curve_container;
    
    if ~any(isnan(cell2mat(curve_container(5,:))))       %Are there test functions? If not - use the multipliers
        tmp     = cell2mat(curve_container(5,:));
        idx     = find(diff(sign(tmp),1,2)); 

        if ~isempty(idx)    %Is there a sign change in the test functions? If not - use the multipliers
            if ~(numel(idx)==1); warning('More than one bifurcation point detected in the current interval. Searching for the first one. Consider lowering the step width and reexamin this parameter interval.'); end
            idx = idx(1);
            stab_fcn = tmp(idx,:);
        else
            stab_fcn = select_by_multipliers(obj,DYN);  %Get stability function for multipliers

        end
    else
            stab_fcn = select_by_multipliers(obj,DYN); %Get stability function for multipliers
    end

end

function stab_fcn = select_by_multipliers(obj,DYN)
    
            multipliers = cell2mat(obj.curve_container(4,:));   %Get the multipliers
            crit_multi = obj.crit_multi(DYN,multipliers);       %Get the critical values (function correctly gives back either real parts of eigenvalues or absolute value of Floquet multipliers -1 or Lyapunov exponents)
            idx     = find(diff(sign(crit_multi),1,2));         %Find sign change
            if ~(numel(idx)==1); warning('Possibly more than one bifurcation point detected in the current interval. Searching for the first one. Consider lowering the step width and reexamin this parameter interval.'); end
            idx = idx(1);
            stab_fcn = crit_multi(idx,:);
    
end

