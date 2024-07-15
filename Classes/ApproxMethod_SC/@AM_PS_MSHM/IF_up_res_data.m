%This methods updates the necessary data of the ApproxMethod object with
%data from the Continuer object

function obj = IF_up_res_data(obj,CON)

    obj.iv = CON.yp(1:(end-1));             %update the current initial condition. Used for the poincare phase condition.

    % for FGM or FD:
    % Check quality of solution and update AM object properties, so that
    % the residual equation can be built up correspondingly. 

end