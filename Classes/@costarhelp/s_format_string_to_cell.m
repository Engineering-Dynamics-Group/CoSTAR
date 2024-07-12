% s_format_string_to_cell formats a string with a given width by breaking it into
% multiple lines.
%
%   @str: The string to format.
%   @width: The maximum width of each line.
%
%   @output: A cell array of strings, with each string being a line of the
%               formatted output.
%
%   @brealtype: Either "characters" or "words" to indicate how the string should be broken

function output = s_format_string_to_cell(str, width, breaktype)


switch breaktype

    case 'character'

        tmp = split(str,'\n'); %Split at the delimiters
        %delete empty rows
        idx = cellfun(@isempty, tmp); %strrim removes leading and trailling white spaces
        % delete empty rows
        tmp(idx) = [];

        % Initialize output cell array
        masteroutput = cell(length(tmp),1);

        for k = 1:length(tmp)

            output = {''};
            strLength = length(tmp{k,1});
            if strLength <= width
                output{end+1,1} = tmp{k,1};
            else
                nLines = ceil(strLength/width);
                 for i = 1:nLines-1
                    startIdx = (i-1)*width + 1;
                    endIdx = i*width;
                    output{end+1,1} = tmp{k,1}(startIdx:endIdx);
                end
                output{end+1,1} = tmp{k,1}((nLines-1)*width+1:end);
            end
            output(1,:) =[];
            masteroutput{k,1} = cellfun(@strtrim,output,'UniformOutput',false);
        end

        output = vertcat(masteroutput{:});

    case 'word'

        tmp = split(str,'\n'); %Split at the delimiters
        %delete empty rows
        idx = cellfun(@isempty, tmp); %strrim removes leading and trailling white spaces
        % delete empty rows
        tmp(idx) = [];


        % Initialize output cell array
        masteroutput = cell(length(tmp),1);
        for k = 1:length(tmp)

            words = strsplit(tmp{k,1}, ' ');
            output = {' '};
            % Loop through words and add to output cell array
            for i = 1:length(words)
                word = words{i};

                % If adding the word exceeds the width, start a new line
                if length(output{end}) + length(word) + 1 > width
                    output{end+1,1} = '';
                end

                % Add the word to the current line
                if isempty(output{end})
                    output{end,1} = word;
                else
                    output{end,1} = [output{end} ' ' word];
                end
            end

            masteroutput{k,1} = cellfun(@strtrim,output,'UniformOutput',false);
        end

        output = vertcat(masteroutput{:});

end

end


