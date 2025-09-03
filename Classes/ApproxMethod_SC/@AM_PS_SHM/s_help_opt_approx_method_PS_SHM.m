%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_text: struct of structs describing the parameters and options for the respective opt_approx_method structure


function help_struct = s_help_opt_approx_method_PS_SHM()


    help_struct.info = 'opt_approx_method --- periodic solution --- Shooting method';

    help_struct.mandatory = [];

    help_struct.optional.solver.value   = '''ode45'', ''ode78'', ''ode89'', \n ''ode23'', ''ode113'', ''ode15s'', \n ''ode23s'', ''ode23t'', ''ode23tb''\n  Default: ''ode45''';
    help_struct.optional.solver.text    = 'Defines the Matlab-proprietary numerical time integration algorithm used for the integration in the shooting step. Use specialised solvers for stiff ODEs (e.g., ode15s).';

    help_struct.optional.n_shoot.value  = 'positive integer \n e.g.: 1, 2, 3, 4, ... \n Default: 2';
    help_struct.optional.n_shoot.text   = 'Number of intervals into which the period of a solution is divided for multiple shooting. n_shoot = 1 is denoted as single shooting method, while n_shoot > 1 represents multiple shooting.';
    
    help_struct.optional.phase_condition.value  = '''poincare'' or ''integral''\n Default: ''poincare''';
    help_struct.optional.phase_condition.text   = ['Defines the phase condition (PC), which is only allowed for autonomous systems. The PC is required to compute the frequency, which is an unknown variable for autonomous systems.\n' ...
                                                   'The Poincare PC is easier to compute, while the Integral PC is said to minimise shifts in the solution profile.'];
    
end