%This function returns the real part of the given eigenvalues. This is a universal function for usage in superclass methods
%
%@obj:              Stability subclass ST_EQ object
%@DYN:              DynamicalSystem object
%@multipliers:      eigenvalues of the current solution point       
%
%@cm:               critical multipliers

function cm = crit_multi(obj,DYN,multipliers)

    cm = real(multipliers);


end