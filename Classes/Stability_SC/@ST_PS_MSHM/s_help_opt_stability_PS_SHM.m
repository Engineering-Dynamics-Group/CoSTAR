%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective opt_approx_method structure


function help_struct = s_help_opt_stability_PS_SHM()


    help_struct.info = 'opt_stability --- periodic --- shooting method';

    help_struct.mandatory = [];
   
    help_struct.optional.iterate_bfp.value       = '''on'' or ''off'' \n Default: ''on''';
    help_struct.optional.iterate_bfp.text        = 'Defines wheather the exact bifurcation point (BFP) is calculated if a BFP is detected.';

    help_struct.optional.max_iter.value          = 'positive integer \n Default: 10';
    help_struct.optional.max_iter.text           = 'Maximum number of iterations to calculate the bifurcation point if iterate_bfp = ''on''.';

    help_struct.optional.abstol_multiplier.value = 'positive scalar \n Default: 1e-4';
    help_struct.optional.abstol_multiplier.text  = 'Absolute tolerance of the decisive multiplier if iterate_bfp = ''on''.';

    help_struct.optional.solver.value = 'string \n e.g.: ''ode45'', ''ode15s'', ... \n Default: ''ode45''';
    help_struct.optional.solver.text = ['Defines the solver to be used for the shooting operation. ' ...
                                        'Every currently available Matlab solver is available for the ' ...
                                        'integration.'];
    

end