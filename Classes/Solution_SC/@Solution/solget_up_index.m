%This function updates option.index and sets it correctly according to
%options. 
%
%
%@obj:          Solution Class object
%@DYN:          DynamicalSystem Class object
%@options:      option structure for the postprocessing of the solution
%class


function options = solget_up_index(obj,DYN,options)


    % If options.mu is set, transform it into solution index
    if isfield(options,'mu')    

        if ischar(options.mu)
            options.index =  1:size(obj.s,2);
        else
            for k = 1:numel(options.mu)
                [~,tmp(k)] = min(abs(obj.mu - options.mu(k)));
            end
            options.index = reshape(unique(tmp),1,[]);

        end

        options = rmfield(options,'mu');        % No longer needed

    end
    

    % If options.index hasn't been set already (by the user or the statement above) or user specified 'all'
    if ~isfield(options,'index') || strcmpi(options.index,'all')

        options.index = 1:size(obj.s,2);      % All solutions are given back

    end


end