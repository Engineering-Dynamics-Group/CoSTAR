%Gatekeeper function for the Continuation  class. In here, all input
%parameters are checked, before further processing. 
%
%@GC:           object of Gatekeeper class 
%@opt_cont:     user supplied option structure for the continuation

function s_CON_gatekeeper(GC,opt_cont)

    opt_cont_mandatory_fieldnames  = {'mu_limit'};      %Mandatory fieldnames in the opt_cont structure
    opt_cont_allowed_fieldnames    = {'step_width','max_cont_step','pred','predictor','subspace','mu_limit','direction','step_control','plot','step_width_limit','step_control_param'}; %Allowed fieldnames in the opt_cont structure

    predictor_allowed_fieldvalues = {'tangent','secant','parable','cubic'};
    subspace_allowed_fieldvalues = {'pseudo-arc','natural','arclength','taxi'};
    step_control_allowed_fieldvalues = {'on','off','corrector_iterations','norm_corrector','angle','combination'};

    %% Check the opt_sol_method structure
  
    GC.check_fields(opt_cont,'opt_cont',opt_cont_mandatory_fieldnames,opt_cont_allowed_fieldnames);
    GC.speak;

    %% system structure: check the entries

    %Check the mandatory fields first (these are definitively present)
    %%%%%%%%%%%%%%%%%%%%
    GC.check_data(opt_cont.mu_limit, 'opt_cont.mu_limit', 'double', [1,2], []);
    GC.speak;
    if opt_cont.mu_limit(1) >= opt_cont.mu_limit(2)     % Check that mu_limit is a vector of increasing values
        GC.error_msg{1,end+1} = append('Your provided limits for the continuation parameter mu are opt_cont.mu_limit = [',num2str(opt_cont.mu_limit(1)),', ',num2str(opt_cont.mu_limit(2)),'].');
        GC.error_msg{1,end+1} = append('However, mu_limit must be a vector of increasing numerical values (mu_limit(2) > mu_limit(1)).');
    end

    %Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%
    if isfield(opt_cont,'step_width');         GC.check_data(opt_cont.step_width,        'opt_cont.step_width',        'double', 'scalar', 'positive'); end
    if isfield(opt_cont,'max_cont_step');      GC.check_data(opt_cont.max_cont_step,     'opt_cont.max_cont_step',     'double', 'scalar', 'positive'); end
    if isfield(opt_cont,'pred');               GC.check_data(opt_cont.pred,              'opt_cont.pred',              'char',   [],       predictor_allowed_fieldvalues); end
    if isfield(opt_cont,'predictor');          GC.check_data(opt_cont.predictor,         'opt_cont.predictor',         'char',   [],       predictor_allowed_fieldvalues); end
    if isfield(opt_cont,'subspace');           GC.check_data(opt_cont.subspace,          'opt_cont.subspace',          'char',   [],       subspace_allowed_fieldvalues); end
    if isfield(opt_cont,'direction');          GC.check_data(opt_cont.direction,         'opt_cont.direction',         'double', 'scalar', [-1,1]); end
    if isfield(opt_cont,'step_control');       GC.check_data(opt_cont.step_control,      'opt_cont.step_control',      'char',   [],       step_control_allowed_fieldvalues); end
    if isfield(opt_cont,'plot');               GC.check_data(opt_cont.plot,              'opt_cont.plot',              'char',   [],       {'on','off'}); end
    if isfield(opt_cont,'step_width_limit');   GC.check_data(opt_cont.step_width_limit,  'opt_cont.step_width_limit',  'double', [1,2],    'positive'); end
    if isfield(opt_cont,'step_control_param'); GC.check_data(opt_cont.step_control_param,'opt_cont.step_control_param','double', [1,2],    'positive'); end

    if isfield(opt_cont,'pred') && isfield(opt_cont,'predictor')
        GC.error_msg{1,end+1} = 'You provided the opt_cont fields ''pred'' and ''predictor''. This is not possible since both fields define the predictor method.';
    end

    GC.speak();

end
