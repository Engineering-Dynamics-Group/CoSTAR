%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_text: struct of structs describing the parameters and options for the respective opt_init structure


function help_struct = s_help_opt_init_QPS_SHM()


    help_struct.info = 'opt_init --- quasiperiodic solution --- Shooting method \n !!! There are two possibilities to define initial values !!! ';

    help_struct.mandatory = [];
    
    help_struct.optional.iv.value            = '[n_char_iv*dim x 1] \n (double) array \n e.g.: [1;2;3;4;5;6;...] ';
    help_struct.optional.iv.text             = ['Possibility 1: Defines the initial value for the nonlinear system solver fsolve (without autonomous frequency(s) and continuation parameter). \n' ...
                                                'If numel(iv)/dim does not match the value of the opt_approx_method field ''n_char'' (or its default value), ''iv'' is interpolated. \n' ...
                                                'dim: state space dimension of the system (without autonomous frequency) \n n_char_iv: number of characteristics of ''iv'''];

    help_struct.optional.ic.value            = '[dim x 1] array \n e.g.: [1;2] \n Default: zeros(dim,1)';
    help_struct.optional.ic.text             = ['Possibility 2: Defines an initial point in state space. Based on this value, a time integration is started to determine appropriate initial values.\n' ...
                                                'If ''ic'' or ''iv'' are not supplied, the default value is used. \n dim: state space dimension of the system (without autonomous frequency).'];
    
    help_struct.optional.tinit.value            = 'scalar value \n e.g.: 2000 \n Default: 10000';
    help_struct.optional.tinit.text             = 'If ''ic'' is supplied, ''tinit'' defines the integration time to reach the attractor.';
    
    help_struct.optional.deltat.value            = 'scalar value \n e.g.: 10000 \n Default: 15000 ';
    help_struct.optional.deltat.text             = ['If ''ic'' is supplied, ''deltat'' defines the integration time on the attractor. \n ' ...
                                                    'Choose a sufficiently large number of deltat to ensure that the manifold will be filled by the trajectory.'];
    
    help_struct.optional.dt.value               = 'scalar value \n e.g.: 0.2 \n Default: 0.1';
    help_struct.optional.dt.text                = ['If ''ic'' is supplied, ''dt'' defines the integration increment on the attractor. \n ' ...
                                                   'Choose a sufficiently small number for ''dt'' to ensure that the manifold will be filled densly by the trajectory.'];

end