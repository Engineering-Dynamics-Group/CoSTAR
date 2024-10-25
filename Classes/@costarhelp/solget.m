% This is the help method for the solget costaropts struct.
% It displays the help message on the command window and gives further options

function solget(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '------------------ CoSTAR Help for solget -------------------\n',...
    '-------------------------------------------------------------\n',...
    '\n', ...
    'The Solution class method solget returns solutions in  \n' ...
    'different solution spaces such as time, hypertime and  \n'...
    'frequency. The corresponding costaropts option structure \n' ...
    'defines the returned solution.'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    if ~isempty(varargin)
        fprintf('You supplied key words to the specify the help call for comtplot options struct. \n However, there is no specification for system');
    end

            fprintf(default_text);
      
            fprintf('\n');
            help_struct = Solution.s_help_solget();             % Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           % Display the help structure

            fprintf('\n\n\n');

end