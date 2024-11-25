% Method getIV generates the initial value for the defined solution method
% from the provided initial condition

function obj = getIV(obj,DYN)

    obj.iv = DYN.opt_init.ic;

end