% This is the Gatekeeper function of the subclass AM_PS_FDM.
% It is a static method of AM_PS_FDM and checks all input parameters before further processing.
% The function is called from the static method s_AM_gatekeeper of the class ApproxMethod.
%
% @GC:                   object of Gatekeeper class
% @system:               user supplied option structure for the system
% @opt_approx_method:    user supplied option structure for the approximation method
% @opt_init:             user supplied option structure for the initial condition

function s_PS_FDM_gatekeeper(GC,system,opt_approx_method,opt_init)


    %% Check the opt_approx_method structure  

    % Check if all mandatory fields are given and check if all given fields are allowed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt_approx_method_mandatory_fieldnames  = {};                                           % Mandatory fieldsnames of the opt_approx_method structure
    opt_approx_method_allowed_fieldnames    = {'n_int','scheme','approx_order','points'};   % Allowed fieldsnames of the opt_approx_method structure
    
    GC.check_fields(opt_approx_method, 'opt_approx_method', opt_approx_method_mandatory_fieldnames, opt_approx_method_allowed_fieldnames);       
    
    GC.speak();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Check the values of all mandatory fields of opt_approx_method first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Check the values of all optional fields of opt_approx_method now 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Field 'n_int'
    if isfield(opt_approx_method,'n_int')
        GC.check_data(opt_approx_method.n_int, 'opt_approx_method.n_int', 'double', 'scalar', []);
        if (mod(opt_approx_method.n_int,1) ~= 0) || (opt_approx_method.n_int < 2)       % Check if the number of intervals is an integer >= 2 (< 2 does not make any sense)
            GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "n_int" is ', num2str(opt_approx_method.n_int), '. However, only positive integers >= 2 are allowed.');
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that n_int >= (approx_order + 1) (if this is not true: at least one point would be taken twice for the approximation, which does not make sense)
        % Here: It can only be checked if 'approx_order' and 'points' are not given. If they are given: The check is done below after their data types and values have been checked
        % if ~isfield(opt_approx_method,'points') && ~isfield(opt_approx_method,'approx_order')
        %     if opt_approx_method.n_int < 7
        %         GC.error_msg{1,end+1} = append('You are using n_int = ', num2str(opt_approx_method.n_int), ' intervals in and an approximation order of 6.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition n_int >= 7 must apply. Please increase n_int or decrease approx_order.');
        %     end
        % end
        % GC.speak();
    end

    % Logical check: Make sure that 'points' is not given by user when defining 'scheme' and/or 'approx_method'
    if isfield(opt_approx_method,'points') && (isfield(opt_approx_method,'approx_order') || isfield(opt_approx_method,'scheme'))
        GC.error_msg{1,end+1} = 'The opt_approx_method structure contains the field "points" in conjunction with "scheme" and/or "approx_order", which is not allowed.';
        GC.speak();
    end

    % Field 'scheme'
    if isfield(opt_approx_method,'scheme')
        GC.check_data(opt_approx_method.scheme, 'opt_approx_method.scheme', 'char', [], {'central','forward','backward'});
        GC.speak();
    end
    
    % Field 'approx_order'
    if isfield(opt_approx_method,'approx_order')
        GC.check_data(opt_approx_method.approx_order, 'opt_approx_method.approx_order', 'double', 'scalar', [])     % It is not checked that the data value is 'positive', because there is a bug in the gatekeeper that turns 0 into a positive number!
        % If 'scheme' is central (either it is the default value (not given by user) or it is explicitly given by user):
        if ~isfield(opt_approx_method,'scheme') || ( isfield(opt_approx_method,'scheme') && strcmpi(opt_approx_method.scheme,'central') )
            if (mod(opt_approx_method.approx_order,2) ~= 0) || (opt_approx_method.approx_order < 2)         % Check if the approximation order is a positive even number
                GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "approx_order" is ', num2str(opt_approx_method.approx_order), '. However, only positive even numbers are allowed when using scheme "central".');
            end
        % In all other cases ('scheme' is either forward or backward and therefore given by user):
        else
            if (mod(opt_approx_method.approx_order,1) ~= 0) || (opt_approx_method.approx_order < 1)         % Check if the approximation order is a positive integer
                GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "approx_order" is ', num2str(opt_approx_method.approx_order), '. However, only positive integers are allowed.');
            end
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that n_int >= (approx_order + 1) (if this is not true: at least one point would be taken twice for the approximation, which does not make sense)
        % if isfield(opt_approx_method,'n_int')                           % If 'n_int' is given by user
        %     if opt_approx_method.n_int < (opt_approx_method.approx_order + 1)
        %         GC.error_msg{1,end+1} = append('You are using n_int = ', num2str(opt_approx_method.n_int), ' intervals and an approximation order of ', num2str(opt_approx_method.approx_order), '.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition (approx_order + 1) <= n_int must apply. Please increase n_int or decrease approx_order.');
        %     end
        % else                                                            % If 'n_int' is not given by user
        %     if 100 < (opt_approx_method.approx_order + 1)
        %         GC.error_msg{1,end+1} = append('You are using n_int = 100 intervals and an approximation order of ', num2str(opt_approx_method.approx_order), '.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition (approx_order + 1) <= 100 must apply. Please increase n_int or decrease approx_order.');
        %     end
        % end
        % GC.speak();
    end

    % Field 'points'
    if isfield(opt_approx_method,'points')
        GC.check_data(opt_approx_method.points, 'opt_approx_method.points', 'double', 'vector', []);
        % Check that all elements of 'points' are unique (the elements must be different pairwise)
        if numel(opt_approx_method.points) ~= numel(unique(opt_approx_method.points))
            GC.error_msg{1,end+1} = 'The field "points" of the opt_approx_method structure contains at least two elements which are equal. However, all elements of "points" must be unique.';
        end
        % Check if the elements of 'points' are integers
        if any(mod(opt_approx_method.points,1))             % any returns true if any element of its input is nonzero -> here: if any element of 'points' is not an integer
            GC.error_msg{1,end+1} = 'The field "points" of the opt_approx_method structure contains at least one non-integer. However, all elements of "points" must be integers.';
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that the number of elements of 'points' is not larger than n_int (n_int equals the number of unknown state space vectors to be solved for)
        % if isfield(opt_approx_method,'n_int')
        %     if length(opt_approx_method.points) > opt_approx_method.n_int
        %          GC.error_msg{1,end+1} = append('The number of elements of the opt_approx_method structure field "points" is ', num2str(length(opt_approx_method.points)), '.');
        %          GC.error_msg{1,end+1} = append('However, the number of elements must not exceed the number of intervals n_int = ', num2str(opt_approx_method.n_int), '.');
        %     end
        % else
        %     if length(opt_approx_method.points) > 100
        %         GC.error_msg{1,end+1} = append('The number of elements of the opt_approx_method structure field "points" is ', num2str(length(opt_approx_method.points)), '.'); 
        %         GC.error_msg{1,end+1} = append('However, the number of elements must not exceed the number of intervals n_int = 100.');
        %     end
        % end
        % GC.speak();
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    %% Check the opt_init structure

    % Check if all mandatory fields are given and check if all given fields are allowed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt_init_mandatory_fieldnames  = {};                                % Mandatory fieldsnames of the opt_init structure
    opt_init_allowed_fieldnames    = {'c0','c1','s1','fdm_sol'};        % Allowed fieldsnames of the opt_init structure

    GC.check_fields(opt_init, 'opt_init', opt_init_mandatory_fieldnames, opt_init_allowed_fieldnames);
    
    GC.speak();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Check the values of all mandatory fields of opt_init first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Check the values of all optional fields of opt_init now 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(opt_init,'c0');  GC.check_data(opt_init.c0, 'opt_init.c0', 'double', [system.dim,1], []); end
    
    if isfield(opt_init,'c1');  GC.check_data(opt_init.c1, 'opt_init.c1', 'double', [system.dim,1], []); end
    
    if isfield(opt_init,'s1');  GC.check_data(opt_init.s1, 'opt_init.s1', 'double', [system.dim,1], []); end

    if isfield(opt_init,'fdm_sol')
        GC.check_data(opt_init.fdm_sol, 'opt_init.fdm_sol', 'double', 'vector', []);
        % Check that fdm_sol is not provided if c0, c1 or s1 is specified
        if isfield(opt_init,'c0') || isfield(opt_init,'c1') || isfield(opt_init,'s1')
            GC.error_msg{1,end+1} = 'The opt_init field "fdm_sol" is not allowed in combination with the opt_init fields "c0", "c1" and "s1".';
        end
        % Check that fdm_sol is a column vector
        if size(opt_init.fdm_sol,1) == 1
            GC.error_msg{1,end+1} = 'The opt_init field "fdm_sol" is a row vector, but should be a column vector.';
        end
        % Check that n_int of fdm_sol (= numel(opt_init.fdm_sol)/system.dim) is an integer
        if mod(numel(opt_init.fdm_sol)/system.dim,1) ~= 0
            GC.error_msg{1,end+1} = append('The size of the opt_init field "fdm_sol" is [', num2str(size(opt_init.fdm_sol,1)), ' x 1]. However, the number of intervals of "fdm_sol" equals');
            GC.error_msg{1,end+1} = append('numel(fdm_sol)/system.dim = ', num2str(numel(opt_init.fdm_sol)/system.dim), ', which is not an integer. This is not allowed.');
        end
        % NOT ANYMORE: The length of fdm_sol must match the discretization, i.e. numel(opt_init.fdm_sol) = (opt_approx_method.n_int)*(system.dim)
        % NOW: If length of fdm_sol does not match discretization, fdm_sol is interpolated. Therefore: Only carry out additional check in lines 159 - 162
        % if isfield(opt_approx_method,'n_int')       % Use the provided value of n_int to check the size of fdm_sol      
        %     if length(opt_init.fdm_sol) ~= system.dim*opt_approx_method.n_int
        %         GC.error_msg{1,end+1} = append('The size of opt_init.fdm_sol is [', num2str(size(opt_init.fdm_sol,1)), ' x 1], but should be [(opt_approx_method.n_int)*(system.dim) x 1] = [', num2str(opt_approx_method.n_int*system.dim), ' x 1].');
        %         GC.error_msg{1,end+1} = 'You provided opt_approx_method.n_int (n_int). If this error occurs, please identify n_int of your provided solution opt_init.fdm_sol and set n_int to its correct value.';
        %     end
        % else                                        % Use the default value of n_int to check the size of fdm_sol                   
        %     if length(opt_init.fdm_sol) ~= system.dim*100
        %         GC.error_msg{1,end+1} = append('The size of opt_init.fdm_sol is [', num2str(size(opt_init.fdm_sol,1)), ' x 1], but should be [(opt_approx_method.n_int)*(system.dim) x 1] = [', num2str(100*system.dim), ' x 1].');
        %         GC.error_msg{1,end+1} = 'You did not provide opt_approx_method.n_int (n_int) and opt_init.fdm_sol does not match the default value of n_int. If this error occurs, please identify n_int of your provided solution opt_init.fdm_sol and set n_int to its correct value.';
        %     end
        % end
    end

    GC.speak();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end