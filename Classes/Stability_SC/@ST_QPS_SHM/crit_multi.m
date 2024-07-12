%This function returns the values of the Lyapunov Exponent. This is a universal function for usage in superclass methods
%
%@obj:              Stability subclass ST_QPS_SHM object
%@DYN:              DynamicalSystem object
%@multipliers:      Lyapunov Exponents of the current solution point    
%
%@cm:               critical multipliers

function cm = crit_multi(obj,DYN,multipliers)

%     if DYN.n_auto  > 0
%         %In the autonomous case, there is always one multiplier, which is 1 or very close to it due to numerical inaccuracies. This one is not of importance for determining the
%         multipliers = abs(multipliers);
%         [~,idx] = sort(abs(multipliers-1),'ascend');                 %Sorting of the absolute (values-1) this way always puts the stab_indicators closest to 1 at the first position.
%         multipliers(idx(1)) = [];                                    %Deleting the first of the sorted stab_indicators always eliminates the one closest to 1, which is either the one associated to the autonomous frequencies or one close to it.
%     end
    
    cm = multipliers;

end