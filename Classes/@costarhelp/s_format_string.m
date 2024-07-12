% FORMAT_STRING Formats a string to fit within a maximum line width by inserting line breaks.
% 
%   @str: the input string to format
%   @max_width: the maximum width (number of characters) of each line
%
%   @formatted_str: the formatted string with line breaks inserted
    

function formatted_str = s_format_string(str, max_width)

    % Initialize variables
    formatted_str = ''; % formatted output string
    current_line = ''; % current line being constructed
    words = strsplit(str); % split input string into words
    
    % Loop over words
    for i = 1:numel(words)
        word = words{i};
    
        % If adding the current word to the current line would exceed the maximum width,
        % insert a line break and start a new line
        if numel([current_line ' ' word]) > max_width
            formatted_str = [formatted_str,current_line,sprintf('\n')];
            current_line = '';
        end
    
        % Add the current word to the current line
        current_line = [current_line ' ' word];
    end
    
    % Add the final line to the output string
    formatted_str = [formatted_str current_line];
    
    % Remove leading and trailing whitespace
    formatted_str = strtrim(formatted_str);

end
