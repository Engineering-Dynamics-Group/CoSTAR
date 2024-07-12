% Gatekeeper function for the Stability Equilibrium subclass
% In here, all input parameters are checked, before further processing
% This static method is called from the static method s_ST_gatekeeper
%
% @GC:             object of Gatekeeper class
% @system:         user supplied option structure for the system 
% @opt_sol:        user supplied option structure for the solution
% @opt_stability:  user supplied option structure for the stability


function s_ST_EQ_gatekeeper(GC,system,opt_sol,opt_stability)


    % Check if all mandatory fields are given and check if all given fields are allowed
    %%%%%%%%%%%%%%%%%%%%
    opt_stability_mandatory_fieldnames  = {};                                                   % Mandatory fildsnames in the opt_stability structure
    opt_stability_allowed_fieldnames    = {'abstol_multiplier','max_iter','iterate_bfp'};       % Allowed fieldsnames in the opt_stability structure
    
    GC.check_fields(opt_stability,'opt_stability',opt_stability_mandatory_fieldnames,opt_stability_allowed_fieldnames);
    
    GC.speak();


    % Check the mandatory fields first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%


    % Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%
    if isfield(opt_stability,'abstol_multiplier');      GC.check_data(opt_stability.abstol_multiplier,     'opt_stability.abstol_multiplier',   'double', 'scalar',     'positive'); end
    if isfield(opt_stability,'max_iter');               GC.check_data(opt_stability.max_iter,              'opt_stability.max_iter',            'double', 'scalar',     'positive'); end
    if isfield(opt_stability,'iterate_bfp');            GC.check_data(opt_stability.iterate_bfp,           'opt_stability.iterate_bfp',         'char',   [],     {'on','off'}); end

    GC.speak();


end
