% Function reshapes vector z of quasi-periodic shooting algorithm to be
% able to integrate with ode-solver for all characteristics in one step 
%
%FcnWrapperODE2 of Quasiperiodic Shooting
%This function is a method of subclass AM_QPS_SHM. 
%Function reshapes vector z of quasi-periodic shooting algorithm to be
%able to integrate with ode-solver for all characteristics in one step 
%
%@obj:  ApproximationMethod subclass object
%@t:    Two dimensional time [t(1,:);t(2,:)]
%@z:    Vector of state-space variables of all characteristics
%       (unperturbated and perturbated for Jacobian)
%@Fcn:  RHS of ODE-system
%@PHI:  Phase values for each characteristic
%
%@f:    value of RHS of ODE-system for all characteristics
%
function f = FcnWrapperODE2(obj,t,z,Fcn,PHI)

dim = obj.n;                                                                % Get dimension of state-space
n_char = obj.n_char;                                                        % Get number of characteristics

x = reshape(z,[dim,n_char]);                                                % Reshape vector of initial conditions
F = Fcn(t+PHI,x);                                                           % Evaluate RHS of ODE-system
f = reshape(F,[n_char*dim,1]);                                              % Reshape function values to be a one dimensional vector

end