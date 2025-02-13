% Function fun_Jac_wrapper_init for Quasi-Periodic-Shooting algorithm
% Wrapper function for initial solution for fsolve for Quasi-Periodic-Shooting algorithm
% This function is a method of subclass AM_QPS_SHM. This function provides jacobian and 
% returns residual function and natural parameterization subspace
% constraint as well as the extended jacobian to be able to use option
% "SpecifiyObjectiveGradient" of fsolve
%
%@obj:  ApproximationMethod subclass object
%@y:    solution curve vector (possibly containing auto frequency)
%@y:    intial solution curve vector (possibly containing auto frequency)
%@DYN:  DynamicalSystem class object
%
%@f:    Residuum vector of quasi-periodic-shooting plus value of subspace
%       constraint
%@J:    Jacobian Matrix of quasi-periodic-shooting plus Jacobian of
%       subspace constraint
%
function [f,J] = fun_Jac_wrapper_init(obj,y,y0)

n = size(y,1);                                                              % Get size of solution curve vector
%% Build function
[F,dF] = obj.res(y);                                                        % Evaluate residual function of quasi-periodic shooting
f = [F;y(end,1)-y0(end,1)];                                                 % Add natural parametrization subspace-constrain to make residual vector

%% Build Jacobian
dmu = sqrt(eps)*(1+abs(y(end,1)));                                          % Calculate differential for bifurcation parameter
dy = [zeros(n-1,1);dmu];                                                    % Define vector of differential for derivative with respect to mu
[deltaF,~] =  obj.res(y+dy);                                                % Calulate residual of Quasi-periodic shooting for perturbed bifurcation parameter
J = [[dF,(deltaF-F)./dmu];[zeros(1,n-1),1]];                                % Calulate Jacobian for initial solution of quasi-periodic shooting

end