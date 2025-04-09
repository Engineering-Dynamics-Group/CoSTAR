%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective opt_approx_method structure


function help_struct = s_help_opt_stability_PS_SHM()

    help_struct.info = 'opt_stability --- periodic solution --- Shooting method';

    help_struct.mandatory = [];

    help_struct.optional.solver.value            = '''ode45'', ''ode78'', ''ode89'', \n ''ode23'', ''ode113'', ''ode15s'', \n ''ode23s'', ''ode23t'', ''ode23tb''\n  Default: ''ode45''';
    help_struct.optional.solver.text             = ['Defines the Matlab-proprietary numerical time integration algorithm for the shooting operation.\n' ...
                                                    'Note: This field is not allowed when the shooting method is already used as approximation method.'];

    help_struct.optional.n_shoot.value           = 'positive integer \n e.g.: 1, 2, 3, 4, ... \n Default: 5';
    help_struct.optional.n_shoot.text            = ['Number of intervals into which the period of a solution is divided for the multiple shooting operation. ' ...
                                                    'n_shoot = 1 is denoted as single shooting, while n_shoot > 1 represents multiple shooting.\n' ...
                                                    'Note: This field is not allowed when the shooting method is already used as approximation method.'];

    help_struct.optional.iterate_bfp.value       = '''on'' or ''off'' \n Default: ''on''';
    help_struct.optional.iterate_bfp.text        = 'Defines whether the exact bifurcation point (BFP) is calculated if a BFP is detected.';

    help_struct.optional.max_iter.value          = 'positive integer \n Default: 10';
    help_struct.optional.max_iter.text           = 'Maximum number of iterations to calculate the bifurcation point if iterate_bfp = ''on''.';

    help_struct.optional.abstol_multiplier.value = 'positive scalar \n Default: 1e-4';
    help_struct.optional.abstol_multiplier.text  = 'Absolute tolerance of the decisive multiplier if iterate_bfp = ''on''.';
    
end