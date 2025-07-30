%This function provides additional constraints to be checked if a
%bifurcation of an equilibrium solution is detected. This is to check
%whether a fold bifurcation occurs or a transcritcal/pitchfork
%bifurcation occurs instead
%
%@obj:          object of Stability subclass ST_QPS
%@idx0:         Index of bifurcation (1 = BF)
%@J:            Jacobian at bifurcation point
%
%@idx:          Index for occurring bifurcation (1 = BF)

function idx = AdditionalConstraints(obj,idx0,J)

% Abstract method, no additional bifurcations to be checked for yet
idx = idx0;                             

end