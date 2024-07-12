%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective opt_approx_method structure


function help_struct = s_help_opt_approx_method_PS_FGM()


    help_struct.info = 'opt_approx_method --- periodic solution --- Fourier-Galerkin method';

    help_struct.mandatory = [];

    help_struct.optional.n_fft.value            = 'power of two \n  e.g., 2^4, 2^5, 2^6, ... \n  Default: 2^6';
    help_struct.optional.n_fft.text             = 'Number of evaluation points along periodic orbit for the Fast Fourier Transformation. Always use powers of two for ideal FFT performance! A higher value leads to more accuracy and higher numerical cost.';
    
    help_struct.optional.phasecond.value        = '''poincare'', ''int_poincare'' \n  Default: ''poincare''';
    help_struct.optional.phasecond.text         = 'Type of phasecondition for autonomous periodic solutions only(!). The autonomous frequency is in general unknown, which necessitates another equation defined by the phasecondition. Poincare condition is recommend for periodic solutions.';

    help_struct.optional.error_control.value    = '''on'', ''off'' \n  Default: ''on''';
    help_struct.optional.error_control.text     = 'The error_control automatically estimates the error made (in frequency space) in the approximation and adapts the number of harmonics or lowers it.';

    help_struct.optional.error_limit.value      = '[1x2] double array \n  Default: [1e-3,0.1]';
    help_struct.optional.error_limit.text       = 'If errror_control is ''on'': error_limit defines the maximal and the minimal error for the error_control: If the estimated error is above the maximum, more harmonics are used, if the error is below minimum, the number of harmonics is reduced.';

    help_struct.optional.ec_iter_max.value      = 'positive integer  \n  e.g., 1,2,3,....  \n  Default: 10';
    help_struct.optional.ec_iter_max.text       = 'If errror_control is ''on'': For adapting the ansatz function to comply with the error_limit, CoSTAR iteratively increases the number of harmonics to meet the tolerance. ec_iter_max defines the maximal number of iterations to do so.';

    help_struct.optional.n_hh_max.value         = 'positive integer \n  e.g., 2,3,4,... \n  Default: Inf';
    help_struct.optional.n_hh_max.text          = 'If errror_control is ''on'': Defines the maximum number of harmonics, which CoSTAR will use to comply with the error_limit.';

end