% This function is a method of the subclass AM_PS_FDM
% It provides the residuum vector to be solved for y by fsolve as well as the corresponding Jacobian matrix.
% In initial solution: Both of them are assigned to "Fcn", which is handed to fsolve. Moreover, the option "SpecifiyObjectiveGradient" is set true.
% That way, fsolve uses the pre-defined Jacobian matrix instead of calculating J by itself, which decreases the computing time significantly.
%
% @obj: Object of AM_PS_FDM
% @y:   Curve point (solution vector to be solved for by fsolve)
% @y0:  Start vector for fsolve
%
% @F:   "Complete" residuum vector function (equation system to be solved for y by fsolve)
% @J:   Jacobian matrix of F
%
function [F,J] = corr_fun_init_FDM(obj,y,y0)

%% Build the residuum vector
[res,J_res] = obj.res(y);                       % Get the residuum of the finite-difference equation system 

F = [res;                                       % The "complete" residuum vector function is residuum of the finite-difference equation system ...
     y(end)-y0(end)];                           % and the natural subspace constraint


%% Build the Jacobian matrix
J = [J_res;                                     % The Jacobian of res has already been calculated by obj.res(y)
     zeros(1,length(y)-1), 1];                  % This is the derivation of the natural subspace contraint with respect to y


end