%This function converts all fields to lower case character, if they are a char.
%@obj:  Gatekeeper class object
%@S:    struct 

function S = fieldvalues_to_lower(obj,S)

    F = fieldnames(S);

    for k = 1:length(F)

        if isa(S.(F{k}),'char')             %If the field value is a char
               S.(F{k}) = lower(S.(F{k}));  %Convert it to lower case
        end

    end

end