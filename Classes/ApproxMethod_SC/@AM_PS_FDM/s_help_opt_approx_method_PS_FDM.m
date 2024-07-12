% This function is a method of the subclass AM_PS_FDM.
% It is a help function file for the costarhelp functionality.
%
% This help file is identified by its name.
%
% @help_struct: struct of structs describing the parameters and options for the respective opt_approx_method structure

function help_struct = s_help_opt_approx_method_PS_FDM()

    help_struct.info = ['opt_approx_method --- periodic solution --- finite difference method \n\n' ...
                        'Nomenclature used in the following: \n' ...
                        '  i            grid point index \n' ...
                        '  DeltaTheta   hyper-time interval between two consecutive grid points \n' ...
                        '  theta_i      hyper-time theta at theta_i = i * DeltaTheta \n' ...
                        '  z            state space vector of the system \n'];

    help_struct.mandatory = [];

    help_struct.optional.n_int.value        = 'positive integer >= 2 \n e.g.: 25, 50, 100, ... \n Default: 100';
    help_struct.optional.n_int.text         = ['Number of hyper-time intervals DeltaTheta into which the hyper-time period 2*pi is divided (2*pi = n_int * DeltaTheta). \n' ...
                                               'The number of grid points used to discretise one period equals n_int + 1.']; % \n' ...
                                               % 'Note: ''n_int'' must be greater than or equal to (approx_order + 1) or numel(points) \n (depending on whether ''approx_order'' or ''points'' is given, see below).'];
    
    help_struct.optional.scheme.value       = '''central'', ''forward'' or \n ''backward'' \n Default: ''central''';
    help_struct.optional.scheme.text        = ['The scheme refers to the discretisation scheme used to approximate the derivation dz(theta_i)/dtheta. \n' ...
                                               '''central'': The grid points used are symmetrically distributed around theta_i. \n' ...
                                               '''forward'': Only grid points at \n theta >= theta_i are used. \n' ...
                                               '''backward'': Only grid points at \n theta <= theta_i are used. \n' ...
                                               'Note: The property ''points'' is not allowed when specifying ''scheme''.'];

    help_struct.optional.approx_order.value = 'positive integer \n e.g.: 1, 2, 3, 4, ... \n Default: 6';
    help_struct.optional.approx_order.text  = ['approx_order describes the order of approximation of dz(theta_i)/dtheta. This means that the remainding error of the approximation is proportional to (DeltaTheta)^(approx_order). \n' ...
                                               'Note: If scheme = ''central'', ''approx_order'' must be an even positive integer. Moreover, the property ''points'' is not allowed when specifying ''approx_order.''']; % ' ...
                                               % 'and ''approx_order'' must not be greater than (n_int - 1).'];

    help_struct.optional.points.value       = 'integer vector \n e.g.: [-2,-1,0,1] \n  Default: [-3,-2,-1,0,1,2,3]';    % \n \n p = numel(points)
    help_struct.optional.points.text        = ['The vector ''points'' specifies the local grid points which are used to approximate dz(theta_i)/dtheta. \n' ...
                                               'E.g. when ''points'' is set to [-2,-1,0,1], dz(theta_i)/dtheta is approximated by the state space vectors z at the grid points (i-2), (i-1), (i) and (i+1). \n' ...
                                               'Note: All elements of ''points'' must be unique. ' ... % and p must not be greater than n_int. 
                                               'Moreover, the properties ''scheme'' and ''approx_order'' are not allowed when specifying ''points''.'];

end