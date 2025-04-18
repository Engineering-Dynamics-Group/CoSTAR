% Gatekeeper function for the Stability Quasi Periodic Shooting subclass
% In here, all input parameters are checked, before further processing
% This static method is called from the static method s_ST_gatekeeper
%
% @GC:             object of Gatekeeper class
% @system:         user supplied option structure for the system 
% @opt_sol:        user supplied option structure for the solution
% @opt_stability:  user supplied option structure for the stability


function s_ST_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_stability)


    % Check if all mandatory fields are given and check if all given fields are allowed
    %%%%%%%%%%%%%%%%%%%%
    opt_stability_mandatory_fieldnames  = {};                                                                           % Mandatory fildsnames in the opt_stability structure
    opt_stability_allowed_fieldnames    = {'abstol_multiplier','max_iter','iterate_bfp','solver','n_char_st','n_map'};  % Allowed fieldsnames in the opt_stability structure
    
    GC.check_fields(opt_stability,'opt_stability',opt_stability_mandatory_fieldnames,opt_stability_allowed_fieldnames);
    
    GC.speak();


    % Check the mandatory fields first (these are definitively present)
    %%%%%%%%%%%%%%%%%%%%


    % Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%
    if isfield(opt_stability,'abstol_multiplier');  GC.check_data(opt_stability.abstol_multiplier,  'opt_stability.abstol_multiplier',  'double',   'scalar',   'positive');   end
    if isfield(opt_stability,'max_iter');           GC.check_data(opt_stability.max_iter,           'opt_stability.max_iter',           'double',   'scalar',   'positive');   end
    if isfield(opt_stability,'iterate_bfp');        GC.check_data(opt_stability.iterate_bfp,        'opt_stability.iterate_bfp',        'char',     [],         {'on','off'}); end
    GC.speak();
    
    if isfield(opt_stability,'solver')
        if strcmpi(opt_sol.approx_method,'shooting') || strcmpi(opt_sol.approx_method,'shm')
            % 'solver' is not allowed when using shooting method since solver is already set in opt_approx_method
            GC.error_msg{1,end+1} = 'You selected the approximation method "shooting" and provided a solver in opt_stability.';
            GC.error_msg{1,end+1} = 'However, this is not allowed since the solver is already defined via opt_approx_method when using the shooting method.';
            GC.speak();
        end
        solver_allowed_fieldvalues = {'ode45','ode78','ode89','ode23','ode113','ode15s','ode23s','ode23t','ode23tb'};
        GC.check_data(opt_stability.solver, 'opt_stability.solver', 'char', [], solver_allowed_fieldvalues); 
        GC.speak();
    end

    if isfield(opt_stability,'n_char_st')   
        GC.check_data(opt_stability.n_char_st, 'opt_stability.n_char_st', 'double', 'scalar', 'positive');
        if mod(opt_stability.n_char_st,1)~=0
            GC.error_msg{1,end+1} = append('The value of the options field opt_stability.n_char_st is ',num2str(opt_stability.n_char_st),', but it must be an integer!');
        end
        if opt_stability.n_char_st == 0
            GC.error_msg{1,end+1} = append('The value of the options field opt_stability.n_char_st is 0, but it must be a positive integer (> 0)!');
        end
        GC.speak();
    end

    if isfield(opt_stability,'n_map')   
        GC.check_data(opt_stability.n_map, 'opt_stability.n_char_st', 'double', 'scalar', 'positive');
        if mod(opt_stability.n_map,1)~=0
            GC.error_msg{1,end+1} = append('The value of the options field opt_stability.n_map is ',num2str(opt_stability.n_map),', but it must be an integer!');
        end
        if opt_stability.n_map == 0
            GC.error_msg{1,end+1} = append('The value of the options field opt_stability.n_map is 0, but it must be a positive integer (> 0)!');
        end
        GC.speak();
    end


end