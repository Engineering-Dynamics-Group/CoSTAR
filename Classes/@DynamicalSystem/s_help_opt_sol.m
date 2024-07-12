%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective opt_sol structure


function help_struct = s_help_opt_sol()

    help_struct.info = 'opt_sol';
    
    help_struct.mandatory.sol_type.value = '''equilibrium'' or ''eq'', \n ''periodic'' or ''ps'' OR \n ''quasiperiodic'' or ''qps''';
    help_struct.mandatory.sol_type.text  = 'Defines the type of stationary solution to be computed. You can either use the whole word (e.g. ''equilibrium'') or its abbreviation (e.g. ''eq'').';
    
    help_struct.mandatory.cont.value = '''on'' or ''off''';
    help_struct.mandatory.cont.text  = 'Defines whether a continuation is carried out. \n ''on'': solution curve is computed. \n ''off'': single solution point is computed.';

    help_struct.mandatory.stability.value = '''on'' or ''off''';
    help_struct.mandatory.stability.text  = ['Defines whether stability of the solution is computed.\n ' ...
                                            '--> EXCEPTION: Not yet available for quasi-periodic solutions approximated by Fourier-Galerkin method or finite difference method.'];
    
    help_struct.optional.approx_method.value  = '''fourier-galerkin'' or ''fgm'', \n ''shooting'' or ''shm'' OR \n ''finite-difference'' or ''fdm''';
    help_struct.optional.approx_method.text   = ['Defines the approximation method used for the computation of periodic or quasi-periodic solutions. You can either use the whole word (e.g. ''fourier-galerkin'') or its abbreviation (e.g. ''fgm'').\n' ...
                                                '--> Equilibrium: No value allowed.'];

    help_struct.optional.act_param.value  = 'positive integer \n <= length(system.param) \n e.g.: 3';
    help_struct.optional.act_param.text   = ['Defines where the continuation parameter is located within the system.param cell array (see field ''param'' of the "system" options structure) ' ...
                                             'by cell indexing. E.g., if system.param = {D,Omega,kappa} and Omega is the continuation parameter, act_param = 2.'];

    help_struct.optional.non_auto_freq.value  = '[1x1] or [1x2] double array or function handle, \n e.g.: @(mu) [mu, pi*mu]';
    help_struct.optional.non_auto_freq.text   = ['Defines the non-autonomous frequenc(ies) (naf) of a (quasi-)periodic solution. \n' ...
                                                 '--> If the naf is dependent on the continuation parameter mu, the property is a function handle.\n' ...
                                                 '--> If there is no dependence, the property is an array. \n' ...
                                                 '--> The sum of non-autonomous and autonomous frequencies (naf + af) must correspond to the solution type,' ...
                                                 '    i.e. naf + af = 1 for periodic and naf + af = 2 for quasi-periodic solutions.\n' ...
                                                 '--> Equilibrium: No value allowed.'];

    help_struct.optional.auto_freq.value  = '[1x1] or [1x2] double array, \n e.g.: [1, pi]';
    help_struct.optional.auto_freq.text   = ['Defines an initial guess for the autonomous frequenc(ies) (af) of a (quasi-)periodic solution. \n' ...
                                                 '--> The sum of non-autonomous and autonomous frequencies (naf + af) must correspond to the solution type,' ...
                                                 '    i.e. naf + af = 1 for periodic and naf + af = 2 for quasi-periodic solutions.\n' ...
                                                 '--> Equilibrium: No value allowed.'];

end