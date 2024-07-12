% This function is a method of the subclass AM_PS_FDM
% It provides the residuum vector to be solved for y by fsolve as well as the corresponding Jacobian matrix.
% In m_continuation: Both of them are assigned to "Fcn", which is handed to fsolve. Moreover, the option "SpecifiyObjectiveGradient" is set true.
% That way, fsolve uses the pre-defined Jacobian matrix instead of calculating J by itself, which decreases the computing time significantly.
%
% @obj:  Object of AM_PS_FDM
% @y:    Curve point (solution vector to be solved for by fsolve)
% @CONT: Continuation class object
%
% @F:    "Complete" residuum vector function (equation system to be solved for y by fsolve)
% @J:    Jacobian matrix of F
%
function [F,J] = corr_fun_FDM(obj,y,CON)

%% Build the residuum vector
[res,J_res] = obj.res(y);                       % Get the residuum of the finite-difference equation system 

F = [res;                                       % The "complete" residuum vector function fun is residuum of the finite-difference equation system ...
     CON.sub_con(y,CON)];                       % and the subspace constraint


%% Build the Jacobian matrix
J = [J_res;                                     % The Jacobian of res has already been calculated by obj.res(y)
     CON.d_sub_con(y,CON)];                     % Get the derivation of the subspace contraint with respect to y


end