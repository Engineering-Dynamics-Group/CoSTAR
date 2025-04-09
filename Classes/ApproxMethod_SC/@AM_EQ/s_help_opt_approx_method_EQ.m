%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_text: struct of structs describing the parameters and options for the respective opt_approx_method structure


function help_struct = s_help_opt_approx_method_EQ()


    help_struct.info = ['opt_approx_method --- equilibrium solution \n\n', ...
                        ' To calculate equilibrium solutions, no approximation method is required.\n Thus, there are no options to set.'];

    help_struct.mandatory = [];
    help_struct.optional = [];
    

end