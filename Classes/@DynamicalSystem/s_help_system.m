%This is the help function file for the costarhelp functionality. 
%
%This help file is identified by its name. 
%
%@help_struct: struct of structs describing the parameters and options for the respective system structure


function help_struct = s_help_system()

    help_struct.info = 'system';

    help_struct.mandatory.order.value = '0 or 1';
    help_struct.mandatory.order.text = ['Defines the order of the ODE whose stationary solution is approximated. ODEs of second order are not implemented yet.\n' ...
                                         '--> Equilibria: order = 0, since equilibria are defined by algebraic equations.\n' ...
                                         '--> (Quasi-)Periodic: order = 1 '];

    help_struct.mandatory.rhs.value = 'function handle \n e.g.: @(t,z,param) ...\n ... duffing(t,z,param)';
    help_struct.mandatory.rhs.text = ['Function handle of the right hand side of the system dz/dt = f(t,z,param). \n' ...
                                       'The cell param contains typically the continuation parameter.\n' ...
                                      '--> Equilibria: Structure of function arguments/variables must be @(z,param).\n' ...
                                      '--> (Quasi-)Periodic: Structure of function arguments/variables must be @(t,z,param).'];
    
    help_struct.mandatory.dim.value = 'positive integer \n e.g.: 1,2,3,...';
    help_struct.mandatory.dim.text = 'Dimension of the system dz/dt = f(t,z,param), e.g. dim = 2 for duffing oscillator.';

    help_struct.optional.param.value = 'cell array \n e.g.: {D,omega0}';
    help_struct.optional.param.text = ['Cell containing the parameters of the system dz/dt = f(t,z,param). This field is mandatory if a solution curve shall be continued, ' ...
                                        'since param contains the active continuation parameter (see field ''act_param'' of the "opt_sol" options structure). ' ...
                                        'When computing single solutions (no continuation), param can be omitted if there are no parameters.'];

    help_struct.optional.info.value = 'string \n e.g. "Periodic solution curve of Duffing system."';
    help_struct.optional.info.text = 'Option for documentation purposes. Further information of the system / computation can be provided here.';


end