%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective contplot option structure


function help_struct = s_help_solget()


    help_struct.info = 'solget';

    help_struct.mandatory.space.value = '''time'', ''hypertime'' \n or ''frequency''';
    help_struct.mandatory.space.text  = ['Defines the solution space in which the solution is returned.\n' ...
                                         '--> time: Solution is returned with respect to time t.\n' ...
                                         '--> hypertime: Solution is returned with respect to hypertime theta (equilibrium: no theta, periodic: theta in [0,2*pi], quasi-periodic: theta in [0,2*pi]^2).\n' ...                  
                                         '--> frequency: Frequency content of a solution  is returned (amplitudes with respect to angular frequency).'];

    help_struct.mandatory.eval.value = '''euclidean'', ''all'' or \n function handle \n e.g.: @(z) z(1:2), \n @(z) z(:,1:2) or\n @(z) z(:,:,1:2)';
    help_struct.mandatory.eval.text  = ['Defines the evaluation of the solution, i.e. the returned quantity. \n' ...
                                         '--> ''euclidean'': The Euclidean norm of all state variables is returned.\n' ...
                                         '--> ''all'': All state variables are returned.\n' ...
                                         '--> function handle: The output of the function handle is returned. The dimension of input argument depends on solution type and solution space:\n' ...
                                         '    --- equilibrium + ''hypertime'': \n 1D array, e.g. @(z) z(1:2).\n' ...
                                         '    --- quasi-periodic + ''hypertime'': \n  3D array, e.g. @(z) z(:,:,1:2).\n' ...
                                         '    --- all other solution types and spaces: \n 2D array, e.g. @(z) z(:,1:2).\n' ...
                                         'The mentioned examples return the first two state variables.'];

    help_struct.optional.index.value = 'scalar, vector (positive \n integers) or ''all'' \n e.g.: 5 or [1:10] \n Default: ''all''';
    help_struct.optional.index.text  = ['Vector of increasing indices specifying the solutions of the continuation that are evaluated.\n' ...
                                        '''all'' evaluates all computed solutions.\n Not allowed in combination with field mu.'];

    help_struct.optional.mu.value = 'scalar, vector (double) \n or ''all'' \n e.g.: 0.5 or [1,2.4] \n Default: ''all''';
    help_struct.optional.mu.text  = ['Defines the values of the continuation parameter mu at which the solutions are evaluated. ' ...
                                     'The solutions that match the given mu-values the closest are taken.\n' ...
                                     'WARNING: This can lead to ambiguity if there are several solutions with the same (or similar) value of mu, which is why the field index is a better choice in general!\n' ...
                                     '''all'' evaluates all computed solutions. \n Not allowed in combination with index.'];

    help_struct.optional.interval.value = '[1x2] array (non-negativ\n double) e.g.: [0,10]\n Default:\n - Periodic: [0,T]\n - Quasi-Periodic: [0,2*pi]';
    help_struct.optional.interval.text  = 'Defines the starting time and end time of the evaluation time interval in which the solution is evaluated.\n Not allowed for ''hypertime'' solution space.\n T: Periodic time';

    help_struct.optional.resolution.value = 'scalar or [1x2] array \n (positive integers) \n e.g.: 100 or [40,50] \n Default: 200 or [200,200]';
    help_struct.optional.resolution.text  = ['Defines the resolution / number of points along each of the time or hypertime axes within the evaluation interval. ' ...
                                             'This determines the discretization of the evaluated solution and affects the accuracy of the output in ''frequency'' solution space, where the resolution must be an even number.\n' ...
                                             'NOTE: [1x2] array only available for quasi-periodic solutions in ''hypertime'' solution space '...
                                             'approximated by finite-differences (to define individual resolutions for each of the two hypertime axis)! Scalars are accepted as well.'];

end