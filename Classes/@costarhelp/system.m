% This is the help method for the option system.

%The file calls the corresponding s_help_xxx static methods from the corresponding classes
%
%@varargin:     optional keywords to further specify the output

function system(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    default_text = [...
    '-------------------------------------------------------------\n',...
    '------------------- CoSTAR Help for system ------------------\n',...
    '-------------------------------------------------------------\n',...
    'The system costaropts struct is used for defining options\n',...
    'specific for the general dynamical system itself.'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
    if ~isempty(varargin)
        fprintf('You supplied key words to the specify the help call for system options struct. \n However, there is no specification for system');
    end

            fprintf(default_text);
      
            fprintf('\n');
            help_struct = DynamicalSystem.s_help_system();          %Get the help struct.
            costarhelp.s_disp_help_text(help_struct);           %Display the help structure

            fprintf('\n\n\n');
    
end
