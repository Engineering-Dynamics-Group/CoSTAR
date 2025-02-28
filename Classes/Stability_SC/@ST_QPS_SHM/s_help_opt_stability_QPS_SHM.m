%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective opt_approx_method structure


function help_struct = s_help_opt_stability_QPS_SHM()

    help_struct.info = 'opt_stability --- quasi-periodic solution --- Shooting method';

    help_struct.mandatory = [];
   
    help_struct.optional.solver.value            = 'string \n e.g.: ''ode45'', ''ode15s'', ... \n Default: ''ode45''';
    help_struct.optional.solver.text             = ['Defines the Matlab-proprietary numerical time integration algorithm for the shooting operation.\n' ...
                                                    'Note: This field is not allowed when the shooting method is already used as approximation method.'];

    help_struct.optional.n_char_st.value         = 'positive integer \n e.g.: 10, 20, 30, ... \n Default: 100';
    help_struct.optional.n_char_st.text          = ['Defines the number of characteristics which are "shot" over the manifold in order to obtain the mapping matrices. '...
                                                    'More characteristics improve the accuracy of the computed Ljapunov exponents. Use a high number for complex manifolds.'];
    
    help_struct.optional.n_map.value             = 'positive integer \n e.g.: 1e3, 1e4, 1e5, ... \n Default: 2e4';
    help_struct.optional.n_map.text              = ['Defines how many times a perturbation is mapped in order to calculate the Ljapunov exponents. '...
                                                    'In total, a perturbation is mapped from t=0 to t=n_map*T1, where T1=2*pi/omega_1. '...
                                                    'A high number improves the accuracy but leads to increased numerical costs.'];

    help_struct.optional.iterate_bfp.value       = '''on'' or ''off'' \n Default: ''on''';
    help_struct.optional.iterate_bfp.text        = 'Defines whether the exact bifurcation point (BFP) is calculated if a BFP is detected.';

    help_struct.optional.max_iter.value          = 'positive integer \n Default: 10';
    help_struct.optional.max_iter.text           = 'Maximum number of iterations to calculate the bifurcation point if iterate_bfp = ''on''.';

    help_struct.optional.abstol_multiplier.value = 'positive scalar \n Default: 1e-4';
    help_struct.optional.abstol_multiplier.text  = 'Absolute tolerance of the decisive multiplier if iterate_bfp = ''on''.';

end