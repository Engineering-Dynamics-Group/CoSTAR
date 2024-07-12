% This is the help method for the option contplot.

%The file calls the corresponding s_help_xxx static methods from the corresponding classes
%
%@varargin:     optional keywords to further specify the output

function contplot(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '------------------ CoSTAR Help for contplot -----------------\n',...
    '-------------------------------------------------------------\n',...
    '\n', ...
    'The Solution class method contplot displays a continuation plot.\n',...
    'The corresponding costaropts option structure defines the plot.'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    if ~isempty(varargin)
        fprintf('You supplied key words to the specify the help call for comtplot options struct. \n However, there is no specification for system');
    end

            fprintf(default_text);
      
            fprintf('\n');
            help_struct = Solution.s_help_contplot();          %Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           %Display the help structure

            fprintf('\n\n\n');

end
