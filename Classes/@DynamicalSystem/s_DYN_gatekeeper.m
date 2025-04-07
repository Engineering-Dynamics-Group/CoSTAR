%Gatekeeper function for the DynamicalSystem class. In here, all input
%parameters are checked, before further processing.
%
%@GC:           object of Gatekeeper class 
%@system:       user supplied structure for the system
%@opt_sol:      user supplied option structure for the solution
%
%ALL FIELD NAMES MUST BE WRITTEN IN SMALL CAPITALS

function s_DYN_gatekeeper(GC,system,opt_sol)

    system_mandatory_fieldnames  = {'order','rhs','dim'};                       %required fieldnames in the options super structure
    system_allowed_fieldnames    = {'order','rhs','dim','param','info'};        %allowed fieldnames in the options super structure

    opt_sol_mandatory_fieldnames  = {'sol_type','cont','stability'};                                                                                    %required fieldnames in the options super structure
    opt_sol_allowed_fieldnames    = {'sol_type','cont','stability','approx_method','act_param','non_auto_freq','auto_freq','display','freq_limit'};     %allowed fieldnames in the options super structure

    sol_type_allowed_fieldvalues   = {'equilibrium','eq','periodic','ps','quasiperiodic','qps'};
    display_allowed_fieldvalues = {'off','final','iter','iter-detailed','step-control','error-control','full'};

    %% Check the system structure
    GC.check_fields(system,'system',system_mandatory_fieldnames,system_allowed_fieldnames);

    %% Check the solution option structure
    GC.check_fields(opt_sol,'opt_sol',opt_sol_mandatory_fieldnames,opt_sol_allowed_fieldnames);

    %% Display error messages up to now
    GC.speak;

    %% system structure: check the entries
    %Check the mandatory fields first (these are definitively present)
    %%%%%%%%%%%%%%%%%%%%
    GC.check_data(system.order,'system.order','double','scalar',[0,1]);
    GC.check_data(system.rhs,'system.rhs','function_handle', [] ,[]);
    GC.check_data(system.dim,'system.dim','double','scalar',[]);
    if mod(system.dim,1)~=0; GC.error_msg{1,end+1} = append('System.dim has the value ',num2str(system.dim),', but must be an integer!');end

    %Check the optional fields now 
    %%%%%%%%%%%%%%%%%%%%
    if isfield(system,'param')
        if isfield(opt_sol,'act_param')         % The field act_param must be supplied if param is supplied (due to param{DYN.act_param} is the residuum functions)
            GC.check_data(system.param,'system.param','cell', {'scalar','vector'},[]); 
        else
            GC.error_msg{1,end+1} = 'You supplied the field system.param, but did not supply an active parameter via opt_sol.act_param. This is not allowed.';
        end
    end
    if isfield(system,'info'); GC.check_data(system.info,'system.info','char',[],[]); end    
    GC.speak;


    %% opt_sol structure: check the entries  
    %Check the mandatory fields first (these are definitively present)
    %%%%%%%%%%%%%%%%%%%%
    GC.check_data(opt_sol.sol_type,'opt_sol.sol_type','char', [] ,sol_type_allowed_fieldvalues);
    GC.check_data(opt_sol.cont,'opt_sol.cont','char',[],{'on','off'});
    GC.check_data(opt_sol.stability,'opt_sol.stability','char',[],{'on','off'});
   
    %Check the optional fields now
    %%%%%%%%%%%%%%%%%%%%
    %Check act_param: If active parameter is supplied, a parameter vector must be supplied, too.
    % THIS LOGICAL DEPENDENCY MUST BE CHECKED HERE, SINCE length(system.param) IS NEEDED
    if isfield(opt_sol,'act_param') 
        if isfield(system,'param')
            GC.check_data(opt_sol.act_param,'opt_sol.act_param','double','scalar',1:length(system.param)); %act_param indicates the position in param vector of the active parameter and must therefore be within the number of elements in param 
        else
            GC.error_msg{1,end+1} = 'You supplied the field opt_sol.act_param, but did not supply a parameter vector via system.param.';
        end
    else                                        % If 'act_param' is not supplied ...
        if strcmpi(opt_sol.cont,'on')           % opt_sol.cont MUST be 'off'
            GC.error_msg{1,end+1} = 'You set opt_sol.cont to "on", which means you want to continue a solution branch. However, you did not supply an active parameter via opt_sol.act_param.';
        end
    end
    % Check display
    if isfield(opt_sol,'display')          
        GC.check_data(opt_sol.display,'opt_sol.display','char',[],display_allowed_fieldvalues); 
    end
    GC.speak; %This statement assures that up to now - everything is cool with the supplied data sets

    
    %% EQUILIBRIUM solution: Check values for that case (this is doubled code to some extend... but maybe a little bit clearer)
    if strcmpi(opt_sol.sol_type,'equilibrium') || strcmpi(opt_sol.sol_type,'eq')

        %Check them again (allowed fields changed... it is easier this way)
        opt_ep_sol_mandatory_fieldnames  = {'sol_type','cont','stability'};
        opt_eq_sol_allowed_fieldnames    = {'sol_type','cont','act_param','stability','display'};                        %allowed fieldsnames in the options super structure
        GC.check_fields(opt_sol,'opt_sol',opt_ep_sol_mandatory_fieldnames,opt_eq_sol_allowed_fieldnames);
        GC.speak();

        %Check values again
        GC.check_data(system.order,'system.order','double', 'scalar' ,0);
        GC.speak();

        %Check number of arguments into the RHS. This ALWAYS has to be 2 since the RHS is called by (z,param) in the residual function
        if nargin(system.rhs) ~= 2
            GC.error_msg{1,end+1} = append(['Your solution type is "',opt_sol.sol_type, '" via opt_sol.sol_type. ',...
                                            'Your right hand side via system.rhs has ',num2str(nargin(system.rhs)),' argument(s), but it needs the arguments (z,param).']); 
        end
        GC.speak();

        %If no approx_method is supplied - the sol_type must be equilibrium
        if ~isfield(opt_sol,'approx_method')
            if ~(strcmpi(opt_sol.sol_type,'equilibrium') || strcmpi(opt_sol.sol_type,'eq'))
                GC.error_msg{1,end+1} = 'If no opt_sol.approx_method is supplied, the opt_sol.sol_type must be "equilibrium" or "eq"';
            end
        end
        GC.speak();

    end


    %% PERIODIC solution: Check values for that case (this is doubled code to some extend... but maybe a little bit clearer)
    if strcmpi(opt_sol.sol_type,'periodic') || strcmpi(opt_sol.sol_type,'ps')

        opt_p_sol_mandatory_fieldnames  = {'sol_type','cont','stability','approx_method'};                                                                  %needed fieldsnames in the options super structure
        opt_p_sol_allowed_fieldnames    = {'sol_type','cont','stability','approx_method','act_param','non_auto_freq','auto_freq','display','freq_limit'};   %allowed fieldsnames in the options super structure

        p_approx_method_allowed_fieldvalues = {'fourier-galerkin','fgm','shooting','shm','finite-difference','fdm'};

        %Check them again (mandatory fields changed... it is easier this way)
        GC.check_fields(opt_sol,'opt_sol',opt_p_sol_mandatory_fieldnames,opt_p_sol_allowed_fieldnames);
        GC.speak();

        %Check values again
        GC.check_data(system.order,'system.order','double', 'scalar' ,1);
        GC.check_data(opt_sol.approx_method,'opt_sol.approx_method','char', [],p_approx_method_allowed_fieldvalues);
        GC.speak();

        %%%%%%%%%%%%%%
        %Check freq_limit
        if isfield(opt_sol,'freq_limit')
            GC.check_data(opt_sol.freq_limit,'opt_sol.freq_limit','double','scalar',[]);
            if opt_sol.freq_limit <= 0
                GC.error_msg{1,end+1} = append('Your provided frequency limit via opt_sol.freq_limit is ',num2str(opt_sol.freq_limit),'. However, only positive (> 0) double values are allowed.');
            end
            GC.speak();
            freq_limit = opt_sol.freq_limit;    % This is used for the frequency check down below
        else
            freq_limit = 1e-4;                  % Default value of the frequency limit (ATTENTION: Also set in DynamicalSystem, but DYN is not available here!)
        end
                
        %%%%%%%%%%%%%%
        %Check non_auto_freq and auto_freq      
        if isfield(opt_sol,'non_auto_freq')
            GC.check_data(opt_sol.non_auto_freq,'opt_sol.non_auto_freq',{'double','function_handle'},'scalar',[]);
            if isfield(system,'param')                  % Gatekeeper has checked above that act_param exists if param exists
                mu = system.param{opt_sol.act_param};   % Get the continuation parameter
            else                                        % No param array exists - this CAN be the case for single solutions (no continuation)
                mu = 0;                                 % mu is meaningless in this case, but it has to be set
            end
            if isa(opt_sol.non_auto_freq,'function_handle')     % If non_auto_freq was given as function handle
                Omega = opt_sol.non_auto_freq(mu);
            else                                                % If non_auto_freq was given as double          
                Omega = opt_sol.non_auto_freq;
            end
            warn_msg = check_freq(freq_limit,Omega);
            if ~isempty(warn_msg)
                GC.error_msg{1,end+1} = 'The initial value of the non-autonomous frequency, which you provided via opt_sol.non_auto_freq(mu),';
                GC.error_msg{1,end+1} = append('equals Omega = ', num2str(Omega), '. This is below the frequency limit of ', num2str(freq_limit,'%.0e'), '.');
            end
        end
        if isfield(opt_sol,'auto_freq')
            GC.check_data(opt_sol.auto_freq,'opt_sol.auto_freq','double','scalar',[]);
            warn_msg = check_freq(freq_limit,opt_sol.auto_freq);
            if ~isempty(warn_msg)
                GC.error_msg{1,end+1} = 'The initial value of the autonomous frequency, which you provided via opt_sol.auto_freq,';
                GC.error_msg{1,end+1} = append('equals omega = ', num2str(opt_sol.auto_freq), '. This is below the frequency limit of ', num2str(freq_limit,'%.0e'), '.');
            end
        end
        GC.speak();

        %Assure that there is either non_auto_freq OR auto_freq 
        if ~xor(isfield(opt_sol,'non_auto_freq'),isfield(opt_sol,'auto_freq')) %true if either non_auto_freq and auto_freq are set or none of the two
            if isfield(opt_sol,'non_auto_freq')&&isfield(opt_sol,'auto_freq') %in combination with xor statement this statement is unambiguous
                GC.error_msg{1,end+1} = 'If the opt_sol.sol_type is "periodic" or "ps", you must either supply  opt_sol.auto_freq OR opt_sol.non_auto_freq. However, you supplied BOTH!';
            else
                GC.error_msg{1,end+1} = 'If the opt_sol.sol_type is "periodic" or "ps", you must either supply  opt_sol.auto_freq OR opt_sol.non_auto_freq. However, you supplied NONE!';
            end
        end
        GC.speak();

        %Assure that enough non_autonomous frequencies and autonomous frequencies are supplied
        n_freq = 0;
        %Get the number of non_autonomous frequencies depending on whether it is a function handle or not
        if isfield(opt_sol,'non_auto_freq')
            if isa(opt_sol.non_auto_freq,'function_handle'); n_freq = numel(opt_sol.non_auto_freq(1)); else; n_freq = numel(opt_sol.non_auto_freq); end     
        end
        %Get the number of autonomous frequencies 
        if isfield(opt_sol,'auto_freq'); n_freq = numel(opt_sol.auto_freq); end
        if ~(n_freq==1)
                GC.error_msg{1,end+1} = append('You supplied via opt_sol.auto_freq or opt_sol.non_auto_freq n = ',num2str(n_freq),' frequncies. However, only one frequency is allowed for periodicity.');
        end
        GC.speak();

        %%%%%%%%%%%%%%
        %Check number of arguments into the RHS. This ALWAYS has to be 3 since the RHS is called by (t,z,param) in the residual function
        if nargin(system.rhs) ~=3
            GC.error_msg{1,end+1} = append(['Your solution type is "',opt_sol.sol_type, '" via opt_sol.sol_type. ',...
                                            'Your right hand side via system.rhs has ',num2str(nargin(system.rhs)),' argument(s), but it needs the arguments (t,z,param).']); 
        end
        GC.speak();

    end


    %% QUASIPERIODIC solution: Check values for that case (this is doubled code to some extend... but maybe a little bit clearer)
    if strcmpi(opt_sol.sol_type,'quasiperiodic') || strcmpi(opt_sol.sol_type,'qps')

        opt_qp_sol_mandatory_fieldnames  = {'sol_type','cont','stability','approx_method'};                                                                 %needed fieldsnames in the options super structure
        opt_qp_sol_allowed_fieldnames    = {'sol_type','cont','stability','approx_method','act_param','non_auto_freq','auto_freq','display','freq_limit'};  %allowed fieldsnames in the options super structure

        qp_approx_method_allowed_fieldvalues = {'fourier-galerkin','fgm','shooting','shm','finite-difference','fdm'};

        %Check them again (mandatory fields changed... it is easier this way)
        GC.check_fields(opt_sol,'opt_sol',opt_qp_sol_mandatory_fieldnames,opt_qp_sol_allowed_fieldnames);   
        GC.speak();

        %Check values again
        GC.check_data(system.order,'system.order','double', 'scalar' ,1);
        GC.check_data(opt_sol.approx_method,'opt_sol.approx_method','char', [],qp_approx_method_allowed_fieldvalues);
        GC.speak();    

        %%%%%%%%%%%%%%
        %Check freq_limit
        if isfield(opt_sol,'freq_limit')
            GC.check_data(opt_sol.freq_limit,'opt_sol.freq_limit','double','scalar',[]);
            if opt_sol.freq_limit <= 0
                GC.error_msg{1,end+1} = append('Your provided frequency limit via opt_sol.freq_limit is ',num2str(opt_sol.freq_limit),'. However, only positive (> 0) double values are allowed.');
            end
            GC.speak();
            freq_limit = opt_sol.freq_limit;    % This is used for the frequency check down below
        else
            freq_limit = 1e-4;                  % Default value of the frequency limit (ATTENTION: Also set in DynamicalSystem, but DYN is not available here!)
        end

        %%%%%%%%%%%%%%%%
        %Check non_auto_freq and auto_freq
        if isfield(opt_sol,'non_auto_freq')
            GC.check_data(opt_sol.non_auto_freq,'opt_sol.non_auto_freq',{'double','function_handle'},{'scalar','vector'},[]);
            if isfield(system,'param')                  % Gatekeeper has checked above that act_param exists if param exists
                mu = system.param{opt_sol.act_param};   % Get the continuation parameter
            else                                        % No param array exists - this CAN be the case for single solutions (no continuation)
                mu = 0;                                 % mu is meaningless in this case, but it has to be set
            end
            if isa(opt_sol.non_auto_freq,'function_handle')     % If non_auto_freq was given as function handle
                Omega = opt_sol.non_auto_freq(mu);
            else                                                % If non_auto_freq was given as double          
                Omega = opt_sol.non_auto_freq;
            end
            warn_msg = check_freq(freq_limit,Omega);
            if ~isempty(warn_msg)
                if isscalar(Omega)
                    GC.error_msg{1,end+1} = 'The initial value of the non-autonomous frequency, which you provided via opt_sol.non_auto_freq(mu),';
                    GC.error_msg{1,end+1} = append('equals Omega = ', num2str(Omega), '. This is below the frequency limit of ', num2str(freq_limit,'%.0e'), '.');
                elseif numel(Omega) == 2
                    GC.error_msg{1,end+1} = 'The initial values of the non-autonomous frequencies, which you provided via opt_sol.non_auto_freq(mu),';
                    GC.error_msg{1,end+1} = append('equal Omega_1 = ', num2str(Omega(1)), ' and Omega_2 = ', num2str(Omega(2)), '. At least one of them is below the frequency limit of ', num2str(freq_limit,'%.0e'), '.');
                end
            end
        end
        if isfield(opt_sol,'auto_freq')
            GC.check_data(opt_sol.auto_freq,'opt_sol.auto_freq','double',{'scalar','vector'},[]);
            warn_msg = check_freq(freq_limit,opt_sol.auto_freq);
            if ~isempty(warn_msg)
                if isscalar(opt_sol.auto_freq)
                    GC.error_msg{1,end+1} = 'The initial value of the autonomous frequency, which you provided via opt_sol.auto_freq,';
                    GC.error_msg{1,end+1} = append('equals omega = ', num2str(opt_sol.auto_freq), '. This is below the frequency limit of ', num2str(freq_limit,'%.0e'), '.');
                elseif numel(opt_sol.auto_freq) == 2
                    GC.error_msg{1,end+1} = 'The initial values of the autonomous frequencies, which you provided via opt_sol.auto_freq,';
                    GC.error_msg{1,end+1} = append('equal omega_1 = ', num2str(opt_sol.auto_freq(1)), ' and omega_2 = ', num2str(opt_sol.auto_freq(2)), '. At least one of them is below the frequency limit of ', num2str(freq_limit,'%.0e'), '.');
                end
            end
        end
        GC.speak();

        %Assure that enough non_autonomous frequencies and autonomous frequencies are supplied
        n_naf = 0;  n_af = 0;
        %Get the number of non_autonomous frequencies depending on whether it is a function handle or not
        if isfield(opt_sol,'non_auto_freq')
            if isa(opt_sol.non_auto_freq,'function_handle'); n_naf = numel(opt_sol.non_auto_freq(1)); else; n_naf = numel(opt_sol.non_auto_freq); end     
        end
        %Get the number of autonomous frequencies 
        if isfield(opt_sol,'auto_freq'); n_af = numel(opt_sol.auto_freq); end
        if ~((n_naf+n_af)==2)
                GC.error_msg{1,end+1} = append('Currently, only 2-D quasi-periodic solution can be computed. However, you supplied via opt_sol.auto_freq or/ and opt_sol.non_auto_freq n = ',num2str(n_naf+n_af),' frequncies.');
        end
        GC.speak();

        %%%%%%%%%%%%%%
        %Check number of arguments into the RHS. This ALWAYS has to be 3 since the RHS is called by (t,z,param) in the residual function
        if nargin(system.rhs) ~=3
            GC.error_msg{1,end+1} = append(['Your solution type is "',opt_sol.sol_type, '" via opt_sol.sol_type. ',...
                                            'Your right hand side via system.rhs has ',num2str(nargin(system.rhs)),' argument(s), but it needs the arguments (t,z,param).']); 
        end
        GC.speak();
        
    end


end