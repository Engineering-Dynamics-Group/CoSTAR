% This function is a method of the subclass AM_PS_FDM.
% It is a help function file for the costarhelp functionality.
%
% This help file is identified by its name.
%
% @help_struct: struct of structs describing the parameters and options for the respective opt_init structure

function help_struct = s_help_opt_init_PS_FDM()

    help_struct.info = ['opt_init --- periodic solution --- finite difference method \n\n' ...
                        'There are two possibilities to set the initial value:\nEither c0, c1 and/or s1 OR fdm_sol can be used. \n\n' ...
                        'Nomenclature used in the following: \n' ...
                        ' dim     dimension of the state space of the system \n' ...
                        ' n_int   number of intervals into which the hyper-time period 2*pi is divided \n' ...
                        ' theta   2*pi-periodic hyper-time ( theta = mod(omega*t, 2*pi) ) \n'...
                        '         ( omega: angular frequency,  t: time )'];

    help_struct.mandatory = [];

    help_struct.optional.c0.value   = '[dim x 1] array \n e.g.: [1;0] \n Default: zeros(dim,1)';
    help_struct.optional.c0.text    = 'Possibility 1: The initial values are calculated using a first order Fourier series. c0 defines the constant Fourier series coefficient.\nThe field fdm_sol is not allowed when providing c0.';

    help_struct.optional.c1.value   = '[dim x 1] array \n e.g.: [1;0] \n Default: zeros(dim,1)';
    help_struct.optional.c1.text    = 'Possibility 1: The initial values are calculated using a first order Fourier series. c1 defines the Fourier series coefficient of the cos(theta) term.\nThe field fdm_sol is not allowed when providing c1.';

    help_struct.optional.s1.value   = '[dim x 1] array \n e.g.: [0;-1] \n Default: zeros(dim,1)';
    help_struct.optional.s1.text    = 'Possibility 1: The initial values are calculated using a first order Fourier series. s1 defines the Fourier series coefficient of the sin(theta) term.\nThe field fdm_sol is not allowed when providing s1.';

    help_struct.optional.fdm_sol.value   = '[dim*n_int_fdm_sol x 1] array \n e.g.: Solution_object.s \n Default: (no default value)';
    help_struct.optional.fdm_sol.text    = ['Possibility 2: The initial value can be taken from an already calculated solution. fdm_sol takes the method solution vector s (stored in Solution_object.s). \n',...
                                            'If the number of intervals n_int, defined in opt_approx_method (or its default value), does not match the number of intervals of fdm_sol (n_int_fdm_sol), the provided solution is interpolated. \n' ...
                                            'Moreover, the fields c0, c1 and s1 are not allowed when using fdm_sol.'];
    % help_struct.optional.fdm_sol.value   = '[dim*n_int_old x 1] array \n e.g.: Solution_object.s \n Default: (no default value) \n \n n_int_old: number of inter-\n vals corresponding to fdm_sol';
    % help_struct.optional.fdm_sol.text    = ['Possibility 2: The initial value can be taken from an already calculated solution. fdm_sol takes the method solution vector s (stored in Solution_object.s). \n',...
    %                                         'If the number of intervals n_int_old, which was used to calculate fdm_sol, is not equal to n_int, fdm_sol is interpolated to match n_int. Moreover, the fields c0, c1 and s1 are not allowed when providing fdm_sol.'];

end