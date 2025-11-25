% This function is a method of the subclass AM_PS_SHM.
% It passes information between the continuation algrithm and the ApproxMethod subclass.
%
% @obj: ApproxMethod subclass AM_PS_SHM object
% @CON: Continuation class object

function obj = IF_up_res_data(obj,var1,DYN)

    if isa(var1,'Continuation')                 % If var1 is an object of Continuation
        obj.iv = var1.yp(1:(end-1));            % Update the current initial condition. Used for the phase condition
    elseif isa(var1,'double')                   % var1 should be a solution vector (type double) in all other cases
        obj.iv = var1;                          % Set iv to given solution vector x0 (only relevant in initial_solution)
    end 

end