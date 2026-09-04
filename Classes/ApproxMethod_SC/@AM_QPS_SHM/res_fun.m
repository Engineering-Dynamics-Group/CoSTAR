% Function fun_Jac_wrapper for Quasi-Periodic-Shooting algorithm
% Wrapper function for fsolve for Quasi-Periodic-Shooting algorithm
% This function is a method of subclass AM_QPS_SHM. This function provides 
% jacobian and returns residual function for chosen subspace constraint as 
% well as the extended jacobian to be able to use option
% "SpecifiyObjectiveGradient" of fsolve
%
%@obj:  ApproximationMethod subclass object
%@y:    solution curve vector (possibly containing auto frequency)
%@y:    intial solution curve vector (possibly containing auto frequency)
%@DYN:  DynamicalSystem class object
%@CONT: Continuation class object
%
%@f:    Residuum vector of quasi-periodic-shooting plus value of subspace
%       constraint
%@J:    Jacobian Matrix of quasi-periodic-shooting plus Jacobian of
%       subspace constraint
%
function [F,J] = res_fun(obj,y,CON)

n = size(y,1);                                                              % Get size of solution curve vector

%% Build function
[g,dg] = obj.res(y);                                                        % Evaluate residual function of quasi-periodic shooting
F = [g;CON.sub_con(y,CON)];                                                 % Add subspace-constrain to make residual vector

%% Build Jacobian
dmu = sqrt(eps)*(1+abs(y(end,1)));                                          % Calculate differential for bifurcation parameter
dy = [zeros(n-1,1);dmu];                                                    % Define vector of differential for derivative with respect to mu
delta_gp =  obj.res(y+dy);                                              % Calulate residual of Quasi-periodic shooting for perturbed bifurcation parameter in positiv direction
delta_gm =  obj.res(y-dy);                                              % Calulate residual of Quasi-periodic shooting for perturbed bifurcation parameter in negativ direction
J = [[dg,(delta_gp-delta_gm)./(2.*dmu)];CON.d_sub_con(y,CON)];              % Calulate Jacobian for initial solution of quasi-periodic shooting

end

