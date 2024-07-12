%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective opt_init structure

function help_struct = s_help_opt_init_QPS_FGM()


     help_struct.info = 'opt_init --- quasiperiodic solution --- Fourier-Galerkin method \n !!! There are two possibilities to define initial values !!! ';


    help_struct.mandatory.c0.value            = '[dim x 1] array \n  e.g., [0;0]';
    help_struct.mandatory.c0.text             = 'First possibility: c0 defines the constant Fourier series coefficient. \n  dim: state space dimension of the system (without autonomous frequency(ies)).';

    help_struct.mandatory.cmatrix.value       = '[dim x n] matrix \n  e.g., [1,2; 3,4]';
    help_struct.mandatory.cmatrix.text        = ['First possibility: cmatrix defines the Fourier series coefficients for cosine terms. \n  dim: state space dimension of the system (without autonomous frequency(ies)). \n  2 <= n <= (n_hh-1): ' ...
                                                 'The number of coefficients in cmatrix can be equal or lower to the number of harmonics (provided in hmatrix) minus 1. \n  However: The ordering must correspond to hmatrix. If less coefficients are provided, ' ...
                                                 'the missing ones are guessed based on an exponential decay of the norm of the coefficients. '];

    help_struct.mandatory.smatrix.value       = '[dim x n] matrix \n  e.g., [5,6; 7,8]';
    help_struct.mandatory.smatrix.text        = ['First possibility: smatrix defines the Fourier series coefficients for sine terms. \n  dim: state space dimension of the system (without autonomous frequency(ies)). \n  2 <= n <= (n_hh-1): ' ...
                                                 'The number of coefficients in cmatrix can be equal or lower to the number of harmonics (provided in hmatrix) minus 1. \n  However: The ordering must correspond to hmatrix. If less coefficients are provided, ' ...
                                                 'the missing ones are guessed based on an exponential decay of the norm of the coefficients. '];

    help_struct.mandatory.fc0.value           = '[dim*(2*n+1) x 1] array  \n  e.g., [0;0;1;2;3;4;5;6;7;8]';
    help_struct.mandatory.fc0.text            = ['Second possibility: fc0 defines the assembled Fourier coefficient vector: fc0 = [c0; cmatrix(:); smatrix(:)]. This is an alternative to defining c0, cmatrix and smatrix and is ideal for starting a continuation from a previously ' ...
                                                 ' found solution space point. \n  dim: state space dimension of the system (without autonomous frequency(ies)). \n  2 <= n <= n_hh: ' ...
                                                 'The number n can be equal or lower to the number of harmonics provided in hmatrix. \n  However: The ordering in fc0 must correspond to hmatrix. If less coefficients are provided, ' ...
                                                 'the missing ones are guessed based on an exponential decay of the norm of the coefficients. '];

    
    help_struct.mandatory.hmatrix.value       = '[2 x n_hh] array \n  e.g., [0,1,0,2,1,3; 0,0,1,1,-1,0].';
    help_struct.mandatory.hmatrix.text        = ['Defines the harmonics to be used in the Fourier series. hmatrix is mandatory for BOTH the first and second possiblity of defining initial values. It must at least contain the harmonics [0,1,0;0,0,1]. The ordering of the harmonics is arbitrary, but must correspond' ...
                                                 'to the orderding of the coefficients in c0,cmatrix, smatrix or fc0 respectively. Obviously, hmatrix must not contain doubles.  \n n_hh: number of harmonics.'];


    help_struct.optional = [];


end