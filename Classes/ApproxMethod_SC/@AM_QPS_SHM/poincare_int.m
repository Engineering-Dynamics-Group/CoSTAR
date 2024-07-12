% Integral poincare phase condition 
% This function is a method of subclass AM_QPS_SHM.
% Calculates the integral poincare phase condition for solutions which 
% require a phase condition (mixed-case, full-autonomous case)
%
%@obj:  ApproximationMethod subclass object
%@F:    Values of solution along characteristics
%@F1:   Values of derivarive of reference solution along characteristics
%@Omega:DynamicalSystem class object
%@Ik:   Integration interval
%
%@P:    Value of integral poincare phase condition
%
function P = poincare_int(obj,F,F1,Omega,Ik)

h0 = sum(F1.*F,3);                                                          % Scalar product for state-space variables
P = trapz(obj.phi(2,:),trapz(Omega(1,1)*Ik,h0,1));                          % Integrate phase condition numerically    

end