%Main function of the gatekeeper class, which checks initially the option
%structure and calls all further static methods
%
%@obj:      Gatekeeper object
%@options:  struct of option structure
%
% The gatekeepers are NOT case sensitiv

function options = m_gatekeeper(obj,options)


    options_mandatory_fieldnames  = {'system','opt_sol','opt_init'};                                                                %mandatory fieldsnames in the options super structure
    options_allowed_fieldnames    = {'system','opt_sol','opt_init','opt_approx_method','opt_cont','opt_stability'};                 %allowed fieldsnames in the options super structure
  
    obj.check_data(options,'options','struct',[],[]);        %Check if options is a struct
    obj.speak; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options        = obj.fieldnames_to_lower(options);                                                      %Set all field names to lower case characters
    obj.check_fields(options,'options',options_mandatory_fieldnames,options_allowed_fieldnames);            %updates the error_msg property of the gatekeeper    
    obj.speak;                                                                                              %Display the error messages saved in GC.error_msg until now 

    %User does not have to supply an opt_cont structure
    if ~isfield(options,'opt_cont'); options.opt_cont = costaropts(); end                                       %options.opt_cont does not have to be supplied
    if ~isfield(options,'opt_approx_method'); options.opt_approx_method = costaropts(); end                     %options.opt_approx_method does not have to be supplied - however... call to AM_Gatekeepers is always done
    if ~isfield(options,'opt_stability'); options.opt_stability = costaropts(); end                             %options.opt_stability does not have to be supplied 
    

    obj.check_data(options.system,'options.system','struct',[],[]);                                         %Check if options.system is a struct
    obj.check_data(options.opt_sol,'options.opt_sol','struct',[],[]);                                       %Check if options.opt_sol is a struct
    obj.check_data(options.opt_init,'options.opt_init','struct',[],[]);                                     %Check if options.opt_init is a struct
    obj.check_data(options.opt_approx_method,'options.opt_approx_method','struct',[],[]);                   %Check if options.opt_approx_method is a struct
    obj.check_data(options.opt_cont,'options.opt_cont','struct',[],[]);                                     %Check if options.opt_cont is a struct
    obj.check_data(options.opt_stability,'options.opt_stability','struct',[],[]);                           %Check if options.opt_stability is a struct

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check if all structs, where generated using the costaropts function
    fn = fieldnames(options);
    for k = 1:length(fn) 
        if isfield(options.(fn{k,1}),'costaropts')
            options.(fn{k,1}) = rmfield(options.(fn{k,1}),'costaropts');     %Not longer needed
        else
            obj.error_msg{1,end+1} = append('The option.',fn{k,1},' structure was not created using the costaropts function, which is mandatory.');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Convert all fieldnames to lower case
    options.system =            obj.fieldnames_to_lower(options.system);                                    %set fieldnames to lower case 
    options.opt_sol =           obj.fieldnames_to_lower(options.opt_sol);                                   %set fieldnames to lower case 
    options.opt_init =          obj.fieldnames_to_lower(options.opt_init);                                  %set fieldnames to lower case 
    options.opt_approx_method = obj.fieldnames_to_lower(options.opt_approx_method);                         %set fieldnames to lower case 
    options.opt_cont =          obj.fieldnames_to_lower(options.opt_cont);                                  %set fieldnames to lower case 
    options.opt_stability =     obj.fieldnames_to_lower(options.opt_stability);                             %set fieldnames to lower case 


    %Convert all string fields to lower case
    options.system =            obj.fieldvalues_to_lower(options.system);                                   %set string fieldvalues to lower case 
    options.opt_sol =           obj.fieldvalues_to_lower(options.opt_sol);                                  %set string fieldvalues  to lower case 
    options.opt_init =          obj.fieldvalues_to_lower(options.opt_init);                                 %set string fieldvalues  to lower case
    options.opt_approx_method = obj.fieldvalues_to_lower(options.opt_approx_method);                        %set string fieldvalues  to lower case 
    options.opt_cont =          obj.fieldvalues_to_lower(options.opt_cont);                                 %set string fieldvalues  to lower case 
    options.opt_stability =     obj.fieldvalues_to_lower(options.opt_stability);                            %set string fieldvalues  to lower case 


    %After these initial checks, we can proceed with the class gatekeepers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call the class-specific gatekeepers

    DynamicalSystem.s_DYN_gatekeeper(obj,options.system,options.opt_sol);                                                  %Checks ONLY the system and opt_sol structure
    ApproxMethod.s_AM_gatekeeper(obj,options.system,options.opt_sol,options.opt_approx_method,options.opt_init);           %Checks ONLY the approx_method and the opt_init structure
    Stability.s_ST_gatekeeper(obj,options.system,options.opt_sol,options.opt_stability);                                   %Checks the opt_stability struct

    if strcmpi(options.opt_sol.cont,'on')
        Continuation.s_CON_gatekeeper(obj,options.opt_cont);                                                %Checks ONLY the cont_method structure %if it is supplied, call the gatekeeper 
    end

    obj.speak;                                                                                              %Display all error messages of the class gatekeepers

end


