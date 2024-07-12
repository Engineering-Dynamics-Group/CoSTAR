% This function is a method of the subclass AM_QPS_FDM.
% It passes information between the continuation algrithm and the ApproxMethod subclass.
%
% @obj: ApproxMethod subclass AM_QPS_FDM object
% @CON: Continuation class object

function obj = IF_up_res_data(obj,CON)  

    obj.iv = CON.yp(1:(end-1));                 % Update the current initial condition. Used for the phase condition

end

