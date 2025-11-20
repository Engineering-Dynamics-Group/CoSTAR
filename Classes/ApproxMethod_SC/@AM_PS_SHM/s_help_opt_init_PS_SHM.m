%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_text: struct of structs describing the parameters and options for the respective opt_init structure


function help_struct = s_help_opt_init_PS_SHM()


    help_struct.info = ['opt_init --- periodic solution --- shooting method \n\n' ...
                        'Nomenclature used in the following: \n' ...
                        '  dim          state space dimension of the system (without autonomous frequency) \n' ...
                        '  n_shoot      number of shooting points \n' ...
                        '  n_shoot_ic   number of shooting points of ''ic'''];

   
    help_struct.mandatory.ic.value = '[dim*n_shoot_ic x 1] array \n  e.g.: [0;1] ';
    help_struct.mandatory.ic.text  = ['Defines the initial condition(s) (i.e. the shooting points) for the shooting method. \n ' ...
                                      '- If numel(ic) = dim (n_shoot_ic = 1), ic is used as initial condition z(t=0) in state space. Missing initial shooting points are obtained via numerical time integration. \n' ...
                                      '- If numel(ic) = dim*n_shoot (n_shoot = n_shoot_ic), ic is directly used as initial guess for the shooting points. \n' ...
                                      '- In all other cases (n_shoot_ic = numel(ic)/dim =/= n_shoot and n_shoot_ic > 1), the state space values stored in ic are interpolated to obtain the desired shooting points. \n' ...
                                      'The initial condition(s) should be as close as possible to point(s) on the periodic orbit.'];
    
    help_struct.optional = [];

end