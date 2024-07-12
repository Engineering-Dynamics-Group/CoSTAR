%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective contplot option structure


function help_struct = s_help_solplot()


    help_struct.info = 'solplot';

    help_struct.mandatory.space.value = '''time'', ''hypertime'', \n ''trajectory'' or ''frequency''';
    help_struct.mandatory.space.text  = ['Defines the solution space in which the solution is to be plotted.\n' ...
                                         '--> ''time'': Solution displayed with respect to time.\n' ...
                                         '--> ''hypertime'': Solution displayed with respect to hypertime (periodic case: 1D manifold, quasi-periodic case: 2D manifold).\n' ...
                                         '--> ''trajectory'': Trajectory of solution displayed in state space or corresponding subspace (2D or 3D respectively). ADDITIONALLY mandatory: Field xaxis for 2D plots, fields xaxis and yaxis for 3D plots (see below).\n' ...                         
                                         '--> ''frequency'': Frequency content of solution is displayed.'];

    help_struct.mandatory.zaxis.value = '''euclidean'', ''all'' or \n function handle \n e.g.: @(z) z(:,1)\n or @(z) z(:,:,1)';
    help_struct.mandatory.zaxis.text  = ['Defines the quantity to be plotted on the ordinate axis (vertical axis). \n' ...
                                         'There are different possibilities depending on chosen solution space:\n' ...
                                         '--> time, hypertime and frequency space: Choose from \n' ...
                                         '     --- ''euclidean'': The Euclidean norm of all state variables is plotted.\n' ...
                                         '     --- ''all'': All state variables are plotted separately.\n' ...
                                         '     --- function handle: Must return a vector valued output for a 2D-array input, e.g.\n @(z) z(:,1) plots the first state variable. ' ...
                                         '                          EXCEPTION for quasi-periodic hypertime plot: Must return a matrix valued output for a 3D-array input, e.g. @(z) z(:,:,1) plots the first state variable.\n' ...
                                         '--> trajectory space: Only function handle allowed. Must return a vector valued output for a 2D-array input, e.g. @(z) z(:,1) plots the first state variable.'];

    help_struct.mandatory.xaxis.value = 'function handle \n e.g.: @(z) z(:,2)';
    help_struct.mandatory.xaxis.text  = ['Defines the quantity to be plotted on the abscissa (horizontal axis). ONLY allowed and mandatory for plots in trajectory solution space (2D and 3D).\n' ...
                                         'The function handle must return a vector valued output for a 2D-array input, e.g.\n @(z) z(:,2) plots the second state variable.'];

    help_struct.mandatory.yaxis.value = 'function handle \n e.g.: @(z) z(:,3)';
    help_struct.mandatory.yaxis.text  = ['Defines the quantity to be plotted on the third axis. ONLY allowed and mandatory for 3D plots in the trajectory solution space.\n' ...
                                         'The function handle must return a vector valued output for a 2D-array input, e.g.\n @(z) z(:,3) plots the third state variable.'];

    help_struct.optional.index.value = 'scalar, vector (positive \n integers) or ''all'' \n e.g.: 5 or [1:10] \n Default: ''all''';
    help_struct.optional.index.text  = 'Vector of increasing indices specifying the solutions of the continuation to be plotted. \n ''all'' plots all computed solutions.';

    help_struct.optional.mu.value = 'scalar, vector (double) \n or ''all'' \n e.g.: 0.5 or [1, 2.4] \n Default: ''all''';
    help_struct.optional.mu.text  = ['Defines the values of the continuation parameter mu at which the solutions are to be plotted.\n' ...
                                     'The solutions that match the given mu-values the closest are displayed.\n' ...
                                     'WARNING: This might lead to ambiguities for overhanging curves. Try using ''index'' instead.\n' ...
                                     '''all'' plots all computed solutions.'];

    help_struct.optional.interval.value = '[1x2] array (non-negativ\n double) e.g.: [0,10]\n Default: [0, 2*pi]';
    help_struct.optional.interval.text  = 'Defines the start and end point of the evaluation time interval. Not allowed for ''hypertime'' solution space.';

    help_struct.optional.resolution.value = 'scalar or [1x2] array \n (positive integers) \n e.g.: 100 or [40 50] \n Default: 200';
    help_struct.optional.resolution.text  = ['Defines the resolution / number of points along each of the time or hypertime axis. ' ...
                                             'This determines the discretization of the plot and affects the accuracy of the frequency plot.\n' ...
                                             'NOTE: [1x2] array (to define individual resolutions for each hypertime axis) only available for hypertime '...
                                             'plots of quasi-periodic solutions approximated by finite-differences (scalars are accepted as well)!'];

    help_struct.optional.figure.value = 'figure handle \n e.g.: gcf \n Default: (no value)';
    help_struct.optional.figure.text  = 'Defines a figure handle of an opened figure (for plotting multiple continuation curves in the same plot).';

    help_struct.optional.color.value = '''r'', ''g'', ''b'', ''c'', ''m'', \n ''y'' or ''k'' \n Default: ''b'' (and further \n colors when plotting multiple\n solutions in one diagram)';
    help_struct.optional.color.text  = 'Defines the color of the  plotted line. Not allowed for ''hypertime'' solution space of quasi-periodic solutions.';

    help_struct.optional.linestyle.value = '''-'', ''--'', '':'' or ''-.'' \n Default: ''-''';
    help_struct.optional.linestyle.text  = 'Defines the linestyle of the  plotted line.';

end