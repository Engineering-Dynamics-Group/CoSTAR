%This is a static method of the costarhelp class and displays a help text
%on the command window. It differentiates between a s_help_xxx file returning
%a single string and a s_help_xxx file returning a help struct.
%
%@myhelp: string or help option struct


function s_disp_help_text(myhelp)


    if isa(myhelp,'struct') %If the help file is a help option struct: Differentiate between mandatory and optional fields and display the options
    
        fprintf('\n');
        fprintf('\n---------------------- CoSTAR help for: ---------------------\n');    
        fprintf(myhelp.info);
        fprintf('\n\n');

        fprintf('\n------------------ Mandatory option fields: -----------------\n');
        fprintf('\n');

        if ~isempty(myhelp.mandatory)
            costarhelp.s_disp_opt_help_struct(myhelp.mandatory);
        else
            fprintf('There are no mandatory option fields \n\n');
        end
    
        fprintf('\n------------------ Optional option fields: ------------------\n');
        fprintf('\n');
    
        if ~isempty(myhelp.optional)  %If it is a simple string: Display the string
            costarhelp.s_disp_opt_help_struct(myhelp.optional);
        else
            fprintf('There are no optional option fields \n');
        end
    
        fprintf('\n-------------------------------------------------------------\n');
    
        
    elseif isa(myhelp,'char')
    
        fprintf('\n');
        fprintf(myhelp);
        fprintf('\n');
        fprintf('-------------------------------------------------------------\n');
      
    end

    
end