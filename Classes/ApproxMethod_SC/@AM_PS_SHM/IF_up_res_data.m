% This function is a method of the subclass AM_PS_SHM.
% It passes information between the continuation algrithm and the ApproxMethod subclass.
%
% @obj: ApproxMethod subclass AM_PS_SHM object
% @CON: Continuation class object

function obj = IF_up_res_data(obj,CON,DYN)

    if DYN.n_auto == 1                          % Stuff for autonomous systems - nothing to do for heteronomous systems

        obj.iv = CON.yp(1:(end-1));             % Save the predictor point -> only needed for Poincare phase condition

    end

end