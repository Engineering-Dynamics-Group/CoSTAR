%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_text: struct of structs describing the parameters and options for the respective opt_init structure


function help_struct = s_help_opt_init_PS_SHM()


    help_struct.info = 'opt_init --- periodic solution --- Shooting method';

   
    help_struct.mandatory.ic.value            = '[dim x 1] array \n  e.g., [0;1] ';
    help_struct.mandatory.ic.text             = ['Defines the initial condition for the shooting method. \n  dim: state space dimension of the system (without autonomous frequency). \n' ...
                                                'Initial condition should be as close as possible to the a point on the periodic orbit.'];
    
    help_struct.optional = [];

end