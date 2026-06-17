% This function is a method of the subclass AM_QPS_FDM.
% It passes information between the continuation algrithm and the ApproxMethod subclass.
%
% @obj:  ApproxMethod subclass AM_QPS_FDM object
% @var1: Continuation class object OR solution vector x whose dimension was updated in the error control

function obj = IF_up_res_data(obj,var,DYN)  

    if isa(var,'Continuation')                  % If var is an object of Continuation
        obj.iv = var.y0(1:(end-1));             % Update the current initial condition. Used for the phase condition
    elseif isa(var,'double')                    % var should be a solution vector (type double) in all other cases
        obj.iv = var(1:(end-1));                % Set iv to given solution vector x0 (only relevant in initial_solution)
    end  

end

