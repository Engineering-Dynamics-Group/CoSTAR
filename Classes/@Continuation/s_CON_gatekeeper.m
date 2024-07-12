%Gatekeeper function for the Continuation  class. In here, all input
%parameters are checked, before further processing. 
%
%@GC:           object of Gatekeeper class 
%@opt_cont:     user supplied option structure for the continuation

function s_CON_gatekeeper(GC,opt_cont)

    opt_cont_mandatory_fieldnames  = {'mu_limit'};                                                                                                                                      %Mandatory fildsnames in the opt_cont structure
    opt_cont_allowed_fieldnames    = {'step_width','max_cont_step','pred','subspace','mu_limit','direction','step_control','plot','display','step_width_limit','step_control_param'};     %Allowed fieldsnames in the opt_cont structure

    pred_allowed_fieldvalues = {'tangent','secant','parable','cubic'};
    subspace_allowed_fieldvalues = {'pseudo-arc','natural','arclength','taxi'};
    step_control_allowed_fieldvalues = {'on','off','corrector_iterations','norm_corrector','angle','combination','pid'};
    display_allowed_fieldvalues = {'on','off','step_control_info'};

    %% Check the opt_sol_method structure
  
    GC.check_fields(opt_cont,'opt_cont',opt_cont_mandatory_fieldnames,opt_cont_allowed_fieldnames);
    GC.speak;

    %% system structure: check the entries

    %Check the mandatory fields first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%
    GC.check_data(opt_cont.mu_limit, 'opt_cont.mu_limit', 'double', [1,2], []);

    %Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%
    if isfield(opt_cont,'step_width');       GC.check_data(opt_cont.step_width,      'opt_cont.step_width',       'double', 'scalar',   'positive'); end
    if isfield(opt_cont,'max_cont_step');    GC.check_data(opt_cont.max_cont_step,   'opt_cont.max_cont_step',    'double', 'scalar',   'positive'); end
    if isfield(opt_cont,'pred');             GC.check_data(opt_cont.pred,            'opt_cont.pred',             'char',   [],         pred_allowed_fieldvalues); end
    if isfield(opt_cont,'subspace');         GC.check_data(opt_cont.subspace,        'opt_cont.subspace',         'char',   [],         subspace_allowed_fieldvalues); end
    if isfield(opt_cont,'direction');        GC.check_data(opt_cont.direction,       'opt_cont.direction',        'double', 'scalar',   [-1,1]); end
    if isfield(opt_cont,'step_control');     GC.check_data(opt_cont.step_control,    'opt_cont.step_control',     'char',   [],         step_control_allowed_fieldvalues); end
    if isfield(opt_cont,'plot');             GC.check_data(opt_cont.plot,            'opt_cont.plot',             'char',   [],         {'on','off'}); end
    if isfield(opt_cont,'display');          GC.check_data(opt_cont.display,         'opt_cont.display',          'char',   [],         display_allowed_fieldvalues); end
    if isfield(opt_cont,'step_width_limit'); GC.check_data(opt_cont.step_width_limit,'opt_cont.step_width_limit', 'double',   [1,2],    'positive'); end
    
    % Check the field 'step_control_param'. This can be either a 1x5 array in case of PID step control or a 1x2 array in all other cases (exception: 'off')
    if isfield(opt_cont,'step_control_param')
        if isfield(opt_cont,'step_control') && strcmpi(opt_cont.step_control,'off')
                % step_control_param is irrelevant if no step control is used and therefore does not have to be checked
        elseif isfield(opt_cont,'step_control') && strcmpi(opt_cont.step_control,'pid')
            GC.check_data(opt_cont.step_control_param, 'opt_cont.step_control_param', 'double', [1,5], 'positive');
        else    % If 'step_control' is given by user, but is not 'pid', or 'step_control' is not given by user (in this case, default method 'angle' is used)
            GC.check_data(opt_cont.step_control_param, 'opt_cont.step_control_param', 'double', [1,2], 'positive');
        end
    end

    GC.speak();

end
