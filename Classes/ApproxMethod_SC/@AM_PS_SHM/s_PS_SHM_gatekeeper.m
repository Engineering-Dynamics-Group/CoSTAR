%Gatekeeper function for the SHM subclass. In here, all input
%parameters are checked, before further processing. This static method is
%called from the static method AM_gatekeeper
%
%@GC:                   object of Gatekeeper class
%@system:               user supplied option structure for the system
%@opt_approx_method:    user supplied option structure for the solution method
%@opt_init:             user supplied option structure for the initial condition

function s_PS_SHM_gatekeeper(GC,system,opt_approx_method,opt_init)

%% Gatekeeper for opt_approx_method
opt_approx_mandatory_fieldnames  = {};                                         %mandatory fieldsnames in the options super structure
opt_approx_allowed_fieldnames    = {'solver'};                                 %allowed fieldsnames in the options super structure

%% Check the opt_approx_method structure
GC.check_fields(opt_approx_method,'opt_approx_method',opt_approx_mandatory_fieldnames,opt_approx_allowed_fieldnames);
GC.speak;

%% system structure: check the entries
solver_allowed_fieldvalues = {'ode45','ode78','ode89','ode23','ode113','ode15s','ode23s','ode23t','ode23tb'};
%Check the mandatory fields first (these are defintively present)
%%%%%%%%%%%%%%%%%%%%

%Check the optional fields now
%%%%%%%%%%%%%%%%%%%%
if isfield(opt_approx_method,'solver'); GC.check_data(opt_approx_method.solver,'opt_approx_method.solver','char',[],solver_allowed_fieldvalues); end

%% Gatekeeper for opt_init
opt_init_mandatory_fieldnames  = {'ic'};                                    %mandatory fieldsnames in the options super structure
opt_init_allowed_fieldnames    = {'ic'};                             %allowed fieldsnames in the options super structure

%% Check the opt_approx_method structure
GC.check_fields(opt_init,'opt_init',opt_init_mandatory_fieldnames,opt_init_allowed_fieldnames);
GC.speak;

%% system structure: check the entries

%Check the mandatory fields first (these are defintively present)
%%%%%%%%%%%%%%%%%%%%
GC.check_data(opt_init.ic,'opt_init.ic','double',[system.dim,1],[]);

end
