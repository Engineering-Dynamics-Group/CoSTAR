% This function is a method of the subclass AM_PS_SHM.
% It passes information between the continuation algrithm and the ApproxMethod subclass.
%
% @obj: ApproxMethod subclass AM_PS_SHM object
% @CON: Continuation class object

function obj = IF_up_res_data(obj,var,DYN)

    if isa(var,'Continuation')                  % If var is an object of Continuation
        obj.iv = var.y0(1:(end-1));             % Update the current initial condition. Used for the phase condition
    elseif isa(var,'double')                    % var should be a solution vector (type double) in all other cases
        obj.iv = var(1:(end-1));                % Set iv to given solution vector x0 (only relevant in initial_solution)
    end

    % If the integral phase condition is used: Get the last computed trajectory (computed in the residuum function)
    if (DYN.n_auto == 1) && strcmpi(obj.phase_condition,'integral')
        obj.Z0 = obj.Z_traj;                    % This is the converged solution and thus the required preceding solution for the phase condition
    end
    
end