%This function converts all fieldnames to lower case character.
%@S: struct 
%It uses the  function S = RenameField(S, Old, New) from fileexchange

function S = fieldnames_to_lower(obj,S)

    Old = fieldnames(S);
    New = lower(Old);

    %% THIS IS NOT MYCODE: Author: Jan Simon, Heidelberg, (C) 2006-2011 (from fileexchange)
    if isempty(S) && isa(S, 'double')  % Accept [] as empty struct without fields
        return;
    end

    Data  = struct2cell(S);
    Field = fieldnames(S);
    if ischar(Old)
        Field(strcmp(Field, Old)) = {New};
    elseif iscellstr(Old)

        for iField = 1:numel(Old)
            Field(strcmp(Field, Old{iField})) = New(iField);
        end
    else
        error(['JSimon:', mfilename, ':BadInputType'], ...
            'Fields must be a string or cell string!');
    end

    S = cell2struct(Data, Field);


end