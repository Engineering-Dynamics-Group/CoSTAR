%Gatekeeper function for the Equilibrium subclass. In here, all input
%parameters are checked, before further processing. This static method is
%called from the static method AM_gatekeeper
%
%@GC:                   object of Gatekeeper class
%@system:               user supplied option structure for the system 
%@opt_approx_method:       user supplied option structure for the solution
%                       method
%@opt_init:             user supplied option structure for the initial
%                       condition

function s_EQ_gatekeeper(GC,system, opt_approx_method,opt_init)

%% Gatekeeper for opt_approx_method
%Nothing to do until now.

%Usual structure:
opt_approxmeth_mandatory_fieldnames  = {};                                    %mandatory fieldsnames in the options super structure
opt_approxmeth_allowed_fieldnames    = {};                                    %allowed fieldsnames in the options super structure

%% Check the opt_approx_method structure
GC.check_fields(opt_approx_method,'opt_approx_method',opt_approxmeth_mandatory_fieldnames,opt_approxmeth_allowed_fieldnames);
GC.speak;

%% system structure: check the entries
%Check the mandatory fields first (these are defintively present)
%%%%%%%%%%%%%%%%%%%%

%Check the optional fields now
%%%%%%%%%%%%%%%%%%%%

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
