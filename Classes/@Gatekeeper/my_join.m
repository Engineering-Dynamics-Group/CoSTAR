%This function is used for the check_data method of the Gatekeeper and
%converts double scalar vectors or 2D matrices of doubles, chars, and cell of chars
%to chars. No cell of cell of chars allowed. 
%
%@obj     object of class Gatekeeper
%@in:     double scalar, vector or 2D matrix, char, cell of char
%@out:    char

function out = my_join(obj,in)

    
    if isa(in,'double')&&numel(size(in))==2 %Only 2D-Matrix supported at the moment
        if length(in) == 1
            out = num2str(in);
        elseif (length(in)>1)&&(numel(in)==length(in))
            out = append('[',strjoin(arrayfun(@num2str, in, 'UniformOutput', 0),', '),']'); 
        else
            out{1,1} = append('[',strjoin(arrayfun(@num2str, in(1,:), 'UniformOutput', 0),', '),';'); 
            for k = 2:(size(in,1)-1)
                out{1,k} = append('',strjoin(arrayfun(@num2str, in(k,:), 'UniformOutput', 0),', '),';'); 
            end
            out{1,end+1} = append('',strjoin(arrayfun(@num2str, in(end,:), 'UniformOutput', 0),', '),']'); 
            out = strjoin(out);
        end

    elseif isa(in,'char')
        out = in;
    elseif isa(in,'cell')
        out = strjoin(in,' , ');
    else
        error('my_join is unsure what to do with this input');
    end



end