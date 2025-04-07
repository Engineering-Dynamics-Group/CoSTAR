%This function returns the values of the Lyapunov Exponent. This is a universal function for usage in superclass methods
%
%@obj:              Stability subclass ST_QPS_SHM object
%@DYN:              DynamicalSystem object
%@multipliers:      Lyapunov Exponents of the current solution point    
%
%@cm:               critical multipliers

function cm = crit_multi(obj,DYN,multipliers)
    
    cm = multipliers;

end