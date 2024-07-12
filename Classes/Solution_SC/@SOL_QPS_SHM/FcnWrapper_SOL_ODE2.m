%Method of SOL_QPS_SHM: This method is used to do the time-integration of
%the quasi-periodic-shooting algorithm in one time-integration
%
%@obj:      obj of class SOL_QPS_SHM < Solution
%@t:        Two-Dimensional time-vector [t(1,:);t(2,:)]
%@z:        Vector of state-space-variables for all characteristics
%@Fcn:      Right-Hand-Side of ODE
%@PHI:      Spacing for characteristics
%
%@f:        value of right-hand side for all characteristics
function f = FcnWrapper_SOL_ODE2(obj,t,z,Fcn,PHI)

dim = obj.n;                                                                % Get dimension of state-space
n_char = obj.n_char;                                                        % Get number of characteristics
x = reshape(z,[dim,n_char]);                                                % reshape vector of initial conditions

F = Fcn(t+PHI,x);                                                           % evaluate function
f = reshape(F,[n_char*dim,1]);                                              % Reshape function values 

end