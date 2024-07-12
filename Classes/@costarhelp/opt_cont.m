%This is the help file for the opt_approx_method structure
%
%The file calls the corresponding s_help_xxx static methods from the corresponding classes
%
%@varargin:     optional keywords to further specify the output

function opt_cont(varargin)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '----------------- CoSTAR Help for opt_cont ------------------\n',...
    '-------------------------------------------------------------\n',...
    'The opt_cont costaropts struct is used for defining options \n',...
    'specific for the continuation of a solution curve. '];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    if ~isempty(varargin)
        fprintf('You supplied key words to the specify the help call for opt_cont. \n However, there is no specification for opt_cont');
    end

            fprintf(default_text);
      
            fprintf('\n');
            help_struct = Continuation.s_help_opt_cont();          %Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           %Display the help structure

            fprintf('\n\n\n');
    
end