%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective contplot option structure


function help_struct = s_help_contplot()

    help_struct.info = 'contplot';
    
    help_struct.mandatory.zaxis.value =  '''min2'', ''mean2'', ''max2'' \n or function handle, e.g.: \n @(z) abs(z(1)) or \n @(z) max(abs(z(:,1))) or \n @(z) max(abs(z(:,:,1)),[],\n ''all'')';
    help_struct.mandatory.zaxis.text  = ['Defines the quantity to be plotted on the ordinate axis of the plot. \n' ...
                                         '--> ''min2'': The minimum of all euclidean vector norms of every point on e.g. a periodic orbit is plotted over the continuation parameter mu. \n' ...
                                         '--> ''mean2'': The mean of all euclidean vector norms of every point on e.g. a periodic orbit is plotted over the continuation parameter mu. \n' ...
                                         '--> ''max2'': The maximum of all euclidean vector norms of every point on e.g. a periodic orbit is plotted over the continuation parameter mu. \n' ...
                                         '--> function handle: Must return a scalar value for a 1D-array (equilibrium solutions), 2D-array (periodic solutions) or 3D-array (quasi-periodic solutions) input. E.g. @(z)abs(z(1)), @(z)max(abs(z(:,1))) or @(z)max(abs(z(:,:,1)),[],''all'') plots the maximum absolute value of the first state variable.'];

    help_struct.optional.index.value = 'scalar, vector (positive \n integer) or ''all'' \n e.g.: 5 or [1:10] \n Default: ''all''';
    help_struct.optional.index.text  = 'Vector of increasing indices specifying the solutions of the continuation to be plotted. \n ''all'' plots all computed solutions.';


    help_struct.optional.resolution.value = 'scalar or [1x2] array \n (positive integers) \n e.g.: 100 or [40,50] \n Default: \n - 200 for periodic solutions\n - 50 or [50,50] for quasi-\n periodic solutions';
    help_struct.optional.resolution.text  = ['Defines the resolution, i.e. the number of points along each of the hypertime axes (periodic: 1 axis, quasi-periodic: 2 axes) ' ...
                                             'which are used for the postprocessing evaluation (see field ''zaxis'' above). This affects the accuracy of the continuation curve except for equilibrium solutions.\n' ...
                                             'For quasi-periodic solutions, a [1x2] array can be specified to define individual resolutions for each axis. However, this is currently only available when using the finite difference method.'];

    help_struct.optional.figure.value = 'figure handle \n e.g.: gcf \n Default: (no value)';
    help_struct.optional.figure.text  = 'Specifies the figure in which to plot the curve. Can be used to plot multiple continuation curves into the same figure.\n When this field is not defined, a new figure is opened.';

    help_struct.optional.color.value = '''r'', ''g'', ''b'', ''c'', ''m'', ''y'',\n ''k'' or [1x3] rgb array \n e.g.: [0,0.5,1] \n Default: ''b'', except \n -> ''r'' for unstable solutions\n -> dark grey when stability\n changes and the bifurcation\n point is not iterated';
    help_struct.optional.color.text  = 'Defines the colour of the  plotted curve.\n When defining a [1x3] array of rgb values, the minimal allowed value for each element is 0 and the maximal allowed value is 1.';

    help_struct.optional.linestyle.value = '''-'', ''--'', '':'' or ''-.'' \n Default: ''-''';
    help_struct.optional.linestyle.text  = 'Defines the line style of the  plotted curve.';

end