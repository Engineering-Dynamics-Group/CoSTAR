% Function updateoptions updates the properties of a class by a given set
% of propteries options = struct('fieldname',value).
%This function is typically called in a constructor of a class.
%Multiple classes call this function.
%
%@obj           Class object, whos properties are to be update
%@opt_struct    Struct to be updated

function obj = updateoptions(obj,opt_struct)

    fields = fieldnames(opt_struct);
    for k=1:length(fields)
        obj.(fields{k,1}) = opt_struct.(fields{k,1});
    end

end