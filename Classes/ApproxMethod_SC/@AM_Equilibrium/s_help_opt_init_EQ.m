%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_text: struct of structs describing the parameters and options for the respective opt_init structure


function help_struct = s_help_opt_init_EQ()


    help_struct.info = 'opt_init --- equilibrium';

   
    help_struct.mandatory.ic.value            = '[dim x 1] array \n  e.g., [0;1] ';
    help_struct.mandatory.ic.text             = ['Defines an intial guess for the first curve point in state space. \n  dim: state space dimension of the system. \n' ...
                                                'Initial condition should be as close as possible to the suspected first curve point.'];
    
    help_struct.optional = [];

end