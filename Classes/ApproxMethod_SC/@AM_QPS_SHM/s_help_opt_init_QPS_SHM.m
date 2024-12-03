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

    help_struct.optional.c0.value   = '[dim x 1] array \n e.g.: [1;0] \n Default: zeros(dim,1)';
    help_struct.optional.c0.text    = 'Possibility 1: The initial values are calculated using a first order multi- dimensional Fourier series. c0 defines the constant Fourier series coefficient.\nThe field fdm_sol is not allowed when providing c0.';

    help_struct.optional.c1_matrix.value   = '[dim x 1], [dim x 2] or \n [dim x 3] array \n e.g.: [1;0], [1,1;0,0] or \n [1,1,1;0,0,0] \n Default: zeros(dim,3)';
    help_struct.optional.c1_matrix.text    = ['Possibility 1: The initial values are calculated using a first order multi- dimensional Fourier series. \n'...
                                              'The first column of c1_matrix defines the coefficient of the cos(theta_1) term. '...
                                              'The second column of c1_matrix defines the coefficient of the cos(theta_2) term. '...
                                              'The third column of c1_matrix defines the coefficient of the cos(theta_1+theta_2) term. '...
                                              'Non-supplied columns are treated as zero vectors.\nThe field fdm_sol is not allowed when providing c1_matrix.'];

    help_struct.optional.s1_matrix.value   = '[dim x 1], [dim x 2] or \n [dim x 3] array \n e.g.: [0;-1], [0,0;-1,-1] or \n [0,0,0;-1,-1,-1] \n Default: zeros(dim,3)';
    help_struct.optional.s1_matrix.text    = ['Possibility 1: The initial values are calculated using a first order multi- dimensional Fourier series. \n'...
                                              'The first column of s1_matrix defines the coefficient of the sin(theta_1) term. '...
                                              'The second column of s1_matrix defines the coefficient of the sin(theta_2) term. '...
                                              'The third column of s1_matrix defines the coefficient of the sin(theta_1+theta_2) term. '...
                                              'Non-supplied columns are treated as zero vectors.\nThe field fdm_sol is not allowed when providing s1_matrix.'];

end