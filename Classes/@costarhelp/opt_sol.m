% This is the help method for the opt_sol costaropts struct.

%The file calls the corresponding s_help_xxx static methods from the corresponding classes
%
%@varargin:     optional keywords to further specify the output

function opt_sol(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '------------------ CoSTAR Help for opt_sol ------------------\n',...
    '-------------------------------------------------------------\n',...
    'The opt_sol costaropts struct is used for defining the options \n',...
    'which determine the general solution and approximation process.'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    if ~isempty(varargin)
        fprintf('You supplied key words to the specify the help call for opt_sol options struct. \n However, there is no specification for opt_sol.');
    end

            fprintf(default_text);
      
            fprintf('\n');
            help_struct = DynamicalSystem.s_help_opt_sol();          %Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           %Display the help structure

            fprintf('\n\n\n');

end