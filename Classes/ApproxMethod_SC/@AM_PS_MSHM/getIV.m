% Method getIV generates the initial value for the defined ApproxMethod method
% from the provided initial condition
function [obj] = getIV(obj,DYN)

if ~isempty(DYN.auto_freq)                                                  %system is (partly) autonomous
    obj.iv = [DYN.opt_init.ic(:);DYN.auto_freq(:)];                         %Append the initial_condition vector, if it is an autonomous system.
else
    obj.iv = DYN.opt_init.ic(:);
end

end