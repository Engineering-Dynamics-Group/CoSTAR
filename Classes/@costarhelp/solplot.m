% This is the help method for the option solplot.

%The file calls the corresponding s_help_xxx static methods from the corresponding classes
%
%@varargin:     optional keywords to further specify the output

function solplot(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '------------------ CoSTAR Help for solplot ------------------\n',...
    '-------------------------------------------------------------\n',...
    '\n', ...
    'The Solution class method solplot displays solutions in \n', ...
    'different solution spaces such as time, hypertime, frequency \n',...
    'or state space. The corresponding costaropts option structure \n', ...
    'defines the plotted solution.'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    if ~isempty(varargin)
        fprintf('You supplied key words to the specify the help call for comtplot options struct. \n However, there is no specification for system');
    end

            fprintf(default_text);
      
            fprintf('\n');
            help_struct = Solution.s_help_solplot();          %Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           %Display the help structure

            fprintf('\n\n\n');
    
end
