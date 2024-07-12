%Gatekeeper function for the AM_template subclass. In here, all input
%parameters are checked, before further processing. This static method is
%called from the static method s_AM_gatekeeper
%
%@GC:                   object of Gatekeeper class 
%@opt_sol_method:       user supplied option structure for the solution
%method
%@opt_init_sol:         user supplied option structure for initializing a
%solution

function s_AM_template_gatekeeper(GC,opt_sol_method,opt_init_sol)

    opt_sol_method_mandatory_fieldnames  = {};                                      %mandatory fieldsnames in the options super structure
    opt_sol_method_allowed_fieldnames    = {};                                      %allowed fieldsnames in the options super structure

    opt_init_sol_mandatory_fieldnames  = {};                                        %mandatory fieldsnames in the options super structure
    opt_init_sol_allowed_fieldnames    = {};                                        %allowed fieldsnames in the options super structure


    %% Check the opt_sol_method structure
  
    GC.check_fields(opt_sol_method,'opt_sol_method',opt_sol_method_mandatory_fieldnames,opt_sol_method_allowed_fieldnames);
    GC.check_fields(opt_init_sol,'opt_init_sol',opt_init_sol_mandatory_fieldnames,opt_init_sol_allowed_fieldnames);
    GC.speak();

    %% check the entries data types
    %Check the mandatory fields first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%

    %Use the GC.check_data function here
    %GC.check_data(opt_sol_method.fieldname,'opt_sol_method.fieldname',allowed data type (e.g. 'double'), allowed data dimension (e.g. [1,2], or 'scalar' or 'vector' or {'scalar','vector'}),allowed data values e.g. [0,1]); end
    % if one of the field cannot be decided, just right '[]'

    %Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%


    %Use the GC.check_data function here

    %% Do logical checks here

end
