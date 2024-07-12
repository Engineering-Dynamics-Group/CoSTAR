% This is the Gatekeeper function of the subclass AM_QPS_FDM.
% It is a static method of AM_QPS_FDM and checks all input parameters before further processing.
% The function is called from the static method s_AM_gatekeeper of the class ApproxMethod.
%
%@GC:                   object of Gatekeeper class
%@system:               user supplied option structure for the system
%@opt_approx_method:    user supplied option structure for the approximation method
%@opt_init:             user supplied option structure for the initial condition

function s_QPS_FDM_gatekeeper(GC,system,opt_approx_method,opt_init)


    %% Check the opt_approx_method structure  

    % Check if all mandatory fields are given and check if all given fields are allowed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt_approx_method_mandatory_fieldnames  = {};                                                   % Mandatory fieldsnames of the opt_approx_method structure
    opt_approx_method_allowed_fieldnames    = {'n_int_1','scheme_1','approx_order_1','points_1',... % Allowed fieldsnames of the opt_approx_method structure
                                               'n_int_2','scheme_2','approx_order_2','points_2',};
    
    GC.check_fields(opt_approx_method, 'opt_approx_method', opt_approx_method_mandatory_fieldnames, opt_approx_method_allowed_fieldnames);       
    
    GC.speak();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Check the values of all mandatory fields of opt_approx_method first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Check the values of all optional fields of opt_approx_method now 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Field 'n_int_1'
    if isfield(opt_approx_method,'n_int_1')
        GC.check_data(opt_approx_method.n_int_1, 'opt_approx_method.n_int_1', 'double', 'scalar', []);
        if (mod(opt_approx_method.n_int_1,1) ~= 0) || (opt_approx_method.n_int_1 < 2)       % Check if the number of intervals is an integer >= 2 (< 2 does not make any sense)
            GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "n_int_1" is ', num2str(opt_approx_method.n_int_1), '. However, only positive integers >= 2 are allowed.');
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that n_int_1 >= (approx_order_1 + 1) (if this is not true: at least one point would be taken twice for the approximation, which does not make sense)
        % Here: It can only be checked if 'approx_order_1' and 'points_1' are not given. If they are given: The check is done below after their data types and values have been checked
        % if ~isfield(opt_approx_method,'approx_order_1') && ~isfield(opt_approx_method,'points_1')
        %     if opt_approx_method.n_int_1 < 7                            
        %         GC.error_msg{1,end+1} = append('You are using n_int_1 = ', num2str(opt_approx_method.n_int_1), ' intervals in and an approximation order of 6 theta_1-direction.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition n_int_1 >= 7 must apply. Please increase n_int_1 or decrease approx_order_1.');
        %     end
        % end
        % GC.speak();
    end

    % Field 'n_int_2'
    if isfield(opt_approx_method,'n_int_2')
        GC.check_data(opt_approx_method.n_int_2, 'opt_approx_method.n_int_2', 'double', 'scalar', []);
        if (mod(opt_approx_method.n_int_2,1) ~= 0) || (opt_approx_method.n_int_2 < 2)       % Check if the number of intervals is an integer >= 2 (< 2 does not make any sense)
            GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "n_int_2" is ', num2str(opt_approx_method.n_int_2), '. However, only positive integers >= 2 are allowed.');
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that n_int_2 >= (approx_order_2 + 1) (if this is not true: at least one point would be taken twice for the approximation, which does not make sense)
        % Here: It can only be checked if 'approx_order_2' and 'points_2' are not given. If they are given: The check is done below after their data types and values have been checked
        % if ~isfield(opt_approx_method,'approx_order_2') && ~isfield(opt_approx_method,'points_2')
        %     if opt_approx_method.n_int_2 < 7                            
        %         GC.error_msg{1,end+1} = append('You are using n_int_2 = ', num2str(opt_approx_method.n_int_2), ' intervals in and an approximation order of 6 theta_2-direction.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition n_int_2 >= 7 must apply. Please increase n_int_2 or decrease approx_order_2.');
        %     end
        % end
        % GC.speak();
    end

    % Logical check: Make sure that 'points_1' ('points_2') is not given by user when defining 'scheme_1' and/or 'approx_method_1' ('scheme_2' and/or 'approx_method_2')
    if isfield(opt_approx_method,'points_1') && (isfield(opt_approx_method,'approx_order_1') || isfield(opt_approx_method,'scheme_1'))
        GC.error_msg{1,end+1} = 'The opt_approx_method structure contains the field "points_1" in conjunction with "scheme_1" and/or "approx_order_1", which is not allowed.';
    end
    if isfield(opt_approx_method,'points_2') && (isfield(opt_approx_method,'approx_order_2') || isfield(opt_approx_method,'scheme_2'))
        GC.error_msg{1,end+1} = 'The opt_approx_method structure contains the field "points_2" in conjunction with "scheme_2" and/or "approx_order_2", which is not allowed.';
    end
    GC.speak();

    % Fields 'scheme_1' and 'scheme_2'
    if isfield(opt_approx_method,'scheme_1')
        GC.check_data(opt_approx_method.scheme_1, 'opt_approx_method.scheme_1', 'char', [], {'central','forward','backward'});
    end
    if isfield(opt_approx_method,'scheme_2')
        GC.check_data(opt_approx_method.scheme_2, 'opt_approx_method.scheme_2', 'char', [], {'central','forward','backward'});   
    end  
    GC.speak();
  
    % Field 'approx_order_1'
    if isfield(opt_approx_method,'approx_order_1')
        GC.check_data(opt_approx_method.approx_order_1, 'opt_approx_method.approx_order_1', 'double', 'scalar', [])     % It is not checked that the data value is 'positive', because there is a bug in the gatekeeper that turns 0 into a positive number!
        % If 'scheme_1' is central (either it is the default value (not given by user) or it is explicitly given by user):
        if ~isfield(opt_approx_method,'scheme_1') || ( isfield(opt_approx_method,'scheme_1') && strcmpi(opt_approx_method.scheme_1,'central') )
            if (mod(opt_approx_method.approx_order_1,2) ~= 0) || (opt_approx_method.approx_order_1 < 2)     % Check if the approximation order is a positive even number
                GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "approx_order_1" is ', num2str(opt_approx_method.approx_order_1), '. However, only positive even numbers are allowed when using scheme "central".');
            end
        % In all other cases ('scheme_1' is either forward or backward and therefore given by user):
        else
            if (mod(opt_approx_method.approx_order_1,1) ~= 0) || (opt_approx_method.approx_order_1 < 1)     % Check if the approximation order is a positive integer
                GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "approx_order_1" is ', num2str(opt_approx_method.approx_order_1), '. However, only positive integers are allowed.');
            end
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that n_int_1 >= (approx_order_1 + 1) (here: approx_order_1 is given by user. The corresponding check if approx_order_1 is not given was already done above)
        % if isfield(opt_approx_method,'n_int_1')                         % If 'n_int_1' is given by user
        %     if opt_approx_method.n_int_1 < (opt_approx_method.approx_order_1 + 1)
        %         GC.error_msg{1,end+1} = append('You are using n_int_1 = ', num2str(opt_approx_method.n_int_1), ' intervals and an approximation order of ', num2str(opt_approx_method.approx_order_1), ' in theta_1-direction.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition (approx_order_1 + 1) <= n_int_1 must apply. Please increase n_int_1 or decrease approx_order_1.');
        %     end
        % else                                                            % If 'n_int_1' is not given by user
        %     if 50 < (opt_approx_method.approx_order_1 + 1)
        %         GC.error_msg{1,end+1} = append('You are using n_int_1 = 50 intervals and an approximation order of ', num2str(opt_approx_method.approx_order_1), ' in theta_1-direction.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition (approx_order_1 + 1) <= 50 must apply. Please increase n_int_1 or decrease approx_order_1.');
        %     end
        % end
        % GC.speak();
    end    

    % Field 'approx_order_2'
    if isfield(opt_approx_method,'approx_order_2')
        GC.check_data(opt_approx_method.approx_order_2, 'opt_approx_method.approx_order_2', 'double', 'scalar', [])     % It is not checked that the data value is 'positive', because there is a bug in the gatekeeper that turns 0 into a positive number!
        % If 'scheme_2' is central (either it is the default value (not given by user) or it is explicitly given by user):
        if ~isfield(opt_approx_method,'scheme_2') || ( isfield(opt_approx_method,'scheme_2') && strcmpi(opt_approx_method.scheme_2,'central') )
            if (mod(opt_approx_method.approx_order_2,2) ~= 0) || (opt_approx_method.approx_order_2 < 2)     % Check if the approximation order is a positive even number
                GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "approx_order_2" is ', num2str(opt_approx_method.approx_order_2), '. However, only positive even numbers are allowed when using scheme "central".');
            end
        % In all other cases ('scheme_2' is either forward or backward and therefore given by user):
        else
            if (mod(opt_approx_method.approx_order_2,1) ~= 0) || (opt_approx_method.approx_order_2 < 1)     % Check if the approximation order is a positive integer
                GC.error_msg{1,end+1} = append('The data value of the opt_approx_method structure field "approx_order_2" is ', num2str(opt_approx_method.approx_order_2), '. However, only positive integers are allowed.');
            end
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that n_int_2 >= (approx_order_2 + 1) (here: approx_order_2 is given by user. The corresponding check if approx_order_2 is not given was already done above)
        % if isfield(opt_approx_method,'n_int_2')                         % If 'n_int_2' is given by user
        %     if opt_approx_method.n_int_2 < (opt_approx_method.approx_order_2 + 1)
        %         GC.error_msg{1,end+1} = append('You are using n_int_2 = ', num2str(opt_approx_method.n_int_2), ' intervals and an approximation order of ', num2str(opt_approx_method.approx_order_2), ' in theta_2-direction.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition (approx_order_2 + 1) <= n_int_2 must apply. Please increase n_int_2 or decrease approx_order_2.');
        %     end
        % else                                                            % If 'n_int_2' is not given by user
        %     if 50 < (opt_approx_method.approx_order_2 + 1)
        %         GC.error_msg{1,end+1} = append('You are using n_int_2 = 50 intervals and an approximation order of ', num2str(opt_approx_method.approx_order_2), ' in theta_2-direction.');
        %         GC.error_msg{1,end+1} = append('This is not allowed because the condition (approx_order_2 + 1) <= 50 must apply. Please increase n_int_2 or decrease approx_order_2.');
        %     end
        % end
        % GC.speak();
    end  

    % Field 'points_1'
    if isfield(opt_approx_method,'points_1')
        GC.check_data(opt_approx_method.points_1, 'opt_approx_method.points_1', 'double', 'vector', []);
        % Check that all elements of 'points_1' are unique (the elements must be different pairwise)
        if numel(opt_approx_method.points_1) ~= numel(unique(opt_approx_method.points_1))
            GC.error_msg{1,end+1} = 'The field "points_1" of the opt_approx_method structure contains at least two elements which are equal. However, all elements of "points_1" must be unique.';
        end
        % Check if the elements of 'points_1' are integers
        if any(mod(opt_approx_method.points_1,1))             % any returns true if any element of its input is nonzero -> here: if any element of 'points_1' is not an integer
            GC.error_msg{1,end+1} = 'The field "points_1" of the opt_approx_method structure contains at least one non-integer. However, all elements of "points_1" must be integers.';
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that the number of elements of 'points_1' is not larger than n_int_1
        % if isfield(opt_approx_method,'n_int_1')
        %     if length(opt_approx_method.points_1) > opt_approx_method.n_int_1
        %          GC.error_msg{1,end+1} = append('The number of elements of the opt_approx_method structure field "points_1" is ', num2str(length(opt_approx_method.points_1)), '.');
        %          GC.error_msg{1,end+1} = append('However, the number of elements must not exceed the number of intervals n_int_1 = ', num2str(opt_approx_method.n_int_1), '.');
        %     end
        % else
        %     if length(opt_approx_method.points_1) > 50
        %         GC.error_msg{1,end+1} = append('The number of elements of the opt_approx_method structure field "points_1" is ', num2str(length(opt_approx_method.points_1)), '.'); 
        %         GC.error_msg{1,end+1} = append('However, the number of elements must not exceed the number of intervals n_int_1 = 50.');
        %     end
        % end
        % GC.speak();
    end

    % Field 'points_2'
    if isfield(opt_approx_method,'points_2')
        GC.check_data(opt_approx_method.points_2, 'opt_approx_method.points_2', 'double', 'vector', []);
        % Check that all elements of 'points_2' are unique (the elements must be different pairwise)
        if numel(opt_approx_method.points_2) ~= numel(unique(opt_approx_method.points_2))
            GC.error_msg{1,end+1} = 'The field "points_2" of the opt_approx_method structure contains at least two elements which are equal. However, all elements of "points_2" must be unique.';
        end
        % Check if the elements of 'points_2' are integers
        if any(mod(opt_approx_method.points_2,1))             % any returns true if any element of its input is nonzero -> here: if any element of 'points_2' is not an integer
            GC.error_msg{1,end+1} = 'The field "points_2" of the opt_approx_method structure contains at least one non-integer. However, all elements of "points_2" must be integers.';
        end
        GC.speak();
        % DEACTIVATED (-> code works nevertheless): Check that the number of elements of 'points_2' is not larger than n_int_2
        % if isfield(opt_approx_method,'n_int_2')
        %     if length(opt_approx_method.points_2) > opt_approx_method.n_int_2
        %          GC.error_msg{1,end+1} = append('The number of elements of the opt_approx_method structure field "points_2" is ', num2str(length(opt_approx_method.points_2)), '.');
        %          GC.error_msg{1,end+1} = append('However, the number of elements must not exceed the number of intervals n_int_2 = ', num2str(opt_approx_method.n_int_2), '.');
        %     end
        % else
        %     if length(opt_approx_method.points_2) > 50
        %         GC.error_msg{1,end+1} = append('The number of elements of the opt_approx_method structure field "points_2" is ', num2str(length(opt_approx_method.points_2)), '.'); 
        %         GC.error_msg{1,end+1} = append('However, the number of elements must not exceed the number of intervals n_int_2 = 50.');
        %     end
        % end
        % GC.speak();
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    %% Check the opt_init structure

    % Check if all mandatory fields are given and check if all given fields are allowed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt_init_mandatory_fieldnames  = {};                                                                                % Mandatory fieldsnames of the opt_init structure
    opt_init_allowed_fieldnames    = {'c0','c1_matrix','s1_matrix','fdm_sol','n_int_1_fdm_sol','n_int_2_fdm_sol'};      % Allowed fieldsnames of the opt_init structure

    GC.check_fields(opt_init, 'opt_init', opt_init_mandatory_fieldnames, opt_init_allowed_fieldnames);
    
    GC.speak();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Check the values of all mandatory fields of opt_init first (these are defintively present)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Check the values of all optional fields of opt_init now 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Field 'c0'
    if isfield(opt_init,'c0')
        GC.check_data(opt_init.c0, 'opt_init.c0', 'double', [system.dim,1], []);
        GC.speak();
    end

    % Field 'c1_matrix'
    if isfield(opt_init,'c1_matrix')  
        GC.check_data(opt_init.c1_matrix, 'opt_init.c1_matrix', 'double', [], []);
        if size(opt_init.c1_matrix,1) ~= system.dim         % Check that the number of rows matches the system dimension
            GC.error_msg{1,end+1} = append('The number of rows of the opt_init structure field "c1_matrix" is ', num2str(size(opt_init.c1_matrix,1)), '.');
            GC.error_msg{1,end+1} = append('However, the number of rows must match the dimension of the system, which is ', num2str(system.dim), '.');
        end
        if size(opt_init.c1_matrix,2) > 3                   % The number of columns can be <= 3 (missing columns are filled by zero vectors)
            GC.error_msg{1,end+1} = append('The number of columns of the opt_init structure field "c1_matrix" is ', num2str(size(opt_init.c1_matrix,2)), '.');
            GC.error_msg{1,end+1} = 'However, the number of columns must not be greater than 3.';
        end
        GC.speak();
    end

    % Field 's1_matrix'
    if isfield(opt_init,'s1_matrix')  
        GC.check_data(opt_init.s1_matrix, 'opt_init.s1_matrix', 'double', [], []);
        if size(opt_init.s1_matrix,1) ~= system.dim         % Check that the number of rows matches the system dimension
            GC.error_msg{1,end+1} = append('The number of rows of the opt_init structure field "s1_matrix" is ', num2str(size(opt_init.s1_matrix,1)), '.');
            GC.error_msg{1,end+1} = append('However, the number of rows must match the dimension of the system, which is ', num2str(system.dim), '.');
        end
        if size(opt_init.s1_matrix,2) > 3                   % The number of columns can be <= 3 (missing columns are filled by zero vectors)
            GC.error_msg{1,end+1} = append('The number of columns of the opt_init structure field "s1_matrix" is ', num2str(size(opt_init.s1_matrix,2)), '.');
            GC.error_msg{1,end+1} = 'However, the number of columns must not be greater than 3.';
        end
        GC.speak();
    end

    % Field 'n_int_1_fdm_sol'
    if isfield(opt_init,'n_int_1_fdm_sol')
        GC.check_data(opt_init.n_int_1_fdm_sol, 'opt_init.n_int_1_fdm_sol', 'double', 'scalar', []);
        if ~isfield(opt_init,'fdm_sol')                                                                 % Check that there is the field 'fdm_sol'
            GC.error_msg{1,end+1} = 'The opt_init field "n_int_1_fdm_sol" is not allowed if you do not provide a solution via the opt_init field "fdm_sol".';
        end
        if isfield(opt_init,'c0') || isfield(opt_init,'c1_matrix') || isfield(opt_init,'s1_matrix')     % Check that 'n_int_1_fdm_sol' is not provided if c0, c1_matrix or s1_matrix is specified
            GC.error_msg{1,end+1} = 'The opt_init field "n_int_1_fdm_sol" is not allowed in combination with the opt_init fields "c0", "c1_matrix" and "s1_matrix".';
        end
        if (mod(opt_init.n_int_1_fdm_sol,1) ~= 0) || (opt_init.n_int_1_fdm_sol < 2)                     % Check if the number of intervals is an integer >= 2 (< 2 does not make any sense)
            GC.error_msg{1,end+1} = append('The data value of the opt_init field "n_int_1_fdm_sol" is ', num2str(opt_init.n_int_1_fdm_sol), '. However, only positive integers >= 2 are allowed.');
        end
        GC.speak();
    end

    % Field 'n_int_2_fdm_sol'
    if isfield(opt_init,'n_int_2_fdm_sol')
        GC.check_data(opt_init.n_int_2_fdm_sol, 'opt_init.n_int_2_fdm_sol', 'double', 'scalar', []);
        if ~isfield(opt_init,'fdm_sol')                                                                 % Check that there is the field 'fdm_sol'
            GC.error_msg{1,end+1} = 'The opt_init field "n_int_2_fdm_sol" is not allowed if you do not provide a solution via the opt_init field "fdm_sol".';
        end
        if isfield(opt_init,'c0') || isfield(opt_init,'c1_matrix') || isfield(opt_init,'s1_matrix')     % Check that 'n_int_2_fdm_sol' is not provided if c0, c1_matrix or s1_matrix is specified
            GC.error_msg{1,end+1} = 'The opt_init field "n_int_2_fdm_sol" is not allowed in combination with the opt_init fields "c0", "c1_matrix" and "s1_matrix".';
        end
        if (mod(opt_init.n_int_2_fdm_sol,1) ~= 0) || (opt_init.n_int_2_fdm_sol < 2)                     % Check if the number of intervals is an integer >= 2 (< 2 does not make any sense)
            GC.error_msg{1,end+1} = append('The data value of the opt_init field "n_int_2_fdm_sol" is ', num2str(opt_init.n_int_2_fdm_sol), '. However, only positive integers >= 2 are allowed.');
        end
        GC.speak();
    end

    % Field 'fdm_sol'
    if isfield(opt_init,'fdm_sol')
        GC.check_data(opt_init.fdm_sol, 'opt_init.fdm_sol', 'double', 'vector', []);
        if isfield(opt_init,'c0') || isfield(opt_init,'c1_matrix') || isfield(opt_init,'s1_matrix')     % Check that fdm_sol is not provided if c0, c1_matrix or s1_matrix is specified
            GC.error_msg{1,end+1} = 'The opt_init field "fdm_sol" is not allowed in combination with the opt_init fields "c0", "c1_matrix" and "s1_matrix".';
        end
        if size(opt_init.fdm_sol,1) == 1                                                                % Check that fdm_sol is a column vector
            GC.error_msg{1,end+1} = 'The opt_init field "fdm_sol" is a row vector, but should be a column vector.';
        end
        GC.speak();
        % Further checks of n_int_1_fdm_sol and n_int_2_fdm_sol depending on which fields are provided: Check if the values match fdm_sol
        if isfield(opt_init,'n_int_1_fdm_sol') && isfield(opt_init,'n_int_2_fdm_sol')       % Check if provided n_int_1_fdm_sol and n_int_2_fdm_sol match the size of fdm_sol      
            if length(opt_init.fdm_sol) ~= (system.dim)*(opt_init.n_int_1_fdm_sol)*(opt_init.n_int_2_fdm_sol)
                GC.error_msg{1,end+1} = append('The size of the opt_init field "fdm_sol" is [', num2str(size(opt_init.fdm_sol,1)), ' x 1]. However, your provided values of the opt_init fields "n_int_1_fdm_sol" and "n_int_2_fdm_sol" do not match "fdm_sol",');
                GC.error_msg{1,end+1} = append('since (system.dim)*(n_int_1_fdm_sol)*(n_int_2_fdm_sol) = ', num2str(opt_init.n_int_1_fdm_sol*opt_init.n_int_2_fdm_sol*system.dim), ' is not equal to numel(fdm_sol) = ', num2str(numel(opt_init.fdm_sol)), '.');
                GC.error_msg{1,end+1} = 'If this error occurs, please identify n_int_1 and/or n_int_2 of your provided solution "fdm_sol" and set their correct values by using the opt_init fields "n_int_1_fdm_sol" and "n_int_2_fdm_sol".';
            end
        elseif isfield(opt_init,'n_int_1_fdm_sol')                                          % Only n_int_1_fdm_sol is provided: Check that n_int_2 of fdm_sol (= numel(opt_init.fdm_sol)/(opt_init.n_int_1_fdm_sol*system.dim) is an integer                
            if mod(numel(opt_init.fdm_sol)/(opt_init.n_int_1_fdm_sol*system.dim),1) ~= 0
                GC.error_msg{1,end+1} = append('The size of the opt_init field "fdm_sol" is [', num2str(size(opt_init.fdm_sol,1)), ' x 1] and you provided the field "n_int_1_fdm_sol". However, the number of intervals n_int_2 of "fdm_sol" equals');
                GC.error_msg{1,end+1} = append('numel(fdm_sol)/(n_int_1_fdm_sol*system.dim) = ', num2str(numel(opt_init.fdm_sol)/(opt_init.n_int_1_fdm_sol*system.dim)),', which is not an integer. This is not allowed.');
                GC.error_msg{1,end+1} = 'If this error occurs, please identify n_int_1 and/or n_int_2 of your provided solution "fdm_sol" and set their correct values by using the opt_init fields "n_int_1_fdm_sol" and "n_int_2_fdm_sol".';
            end
        elseif isfield(opt_init,'n_int_2_fdm_sol')                                          % Only n_int_2_fdm_sol is provided: Check that n_int_1 of fdm_sol (= numel(opt_init.fdm_sol)/(opt_init.n_int_2_fdm_sol*system.dim) is an integer                  
            if mod(numel(opt_init.fdm_sol)/(opt_init.n_int_2_fdm_sol*system.dim),1) ~= 0
                GC.error_msg{1,end+1} = append('The size of the opt_init field "fdm_sol" is [', num2str(size(opt_init.fdm_sol,1)), ' x 1] and you provided the field "n_int_2_fdm_sol". However, the number of intervals n_int_1 of "fdm_sol" equals');
                GC.error_msg{1,end+1} = append('numel(fdm_sol)/(n_int_2_fdm_sol*system.dim) = ', num2str(numel(opt_init.fdm_sol)/(opt_init.n_int_2_fdm_sol*system.dim)),', which is not an integer. This is not allowed.');
                GC.error_msg{1,end+1} = 'If this error occurs, please identify n_int_1 and/or n_int_2 of your provided solution "fdm_sol" and set their correct values by using the opt_init fields "n_int_1_fdm_sol" and "n_int_2_fdm_sol".';
            end
        else                                                                                % Neither n_int_1_fdm_sol nor n_int_2_fdm_sol are provided: Use n_int_1 and n_int_2 of opt_approx_method to check size of fdm_sol
            if isfield(opt_approx_method,'n_int_1') && isfield(opt_approx_method,'n_int_2')     % If n_int_1 and n_int_2 are provided
                n_int_1 = opt_approx_method.n_int_1;    n_int_2 = opt_approx_method.n_int_2;
            elseif isfield(opt_approx_method,'n_int_1')                                         % If only n_int_1 is provided (use default value for n_int_2)
                n_int_1 = opt_approx_method.n_int_1;    n_int_2 = 50;                       
            elseif isfield(opt_approx_method,'n_int_2')                                         % If only n_int_2 is provided (use default value for n_int_1)
                n_int_1 = 50;                           n_int_2 = opt_approx_method.n_int_2;
            else                                                                                % If neither n_int_1 nor n_int_2 is provided (use default values)
                n_int_1 = 50;                           n_int_2 = 50;
            end
            if numel(opt_init.fdm_sol) ~= system.dim*n_int_1*n_int_2
                GC.error_msg{1,end+1} = append('The size of the opt_init field "fdm_sol" is [', num2str(size(opt_init.fdm_sol,1)), ' x 1] and you provided no number of intervals belonging to "fdm_sol" (the fields "n_int_1_fdm_sol" and/or "n_int_2_fdm_sol").');
                GC.error_msg{1,end+1} = 'However, the values of the opt_approx_method fields "n_int_1" and "n_int_2" (either your provided values or their default values) do not match "fdm_sol",';
                GC.error_msg{1,end+1} = append('since (system.dim)*(n_int_1)*(n_int_2) = ', num2str(system.dim*n_int_1*n_int_2), ' is not equal to numel(fdm_sol) = ', num2str(numel(opt_init.fdm_sol)), '.');
                GC.error_msg{1,end+1} = 'If this error occurs, please identify n_int_1 and/or n_int_2 of your provided solution "fdm_sol" and set their correct values by using the opt_init fields "n_int_1_fdm_sol" and "n_int_2_fdm_sol".';
            end
        end
        % NOT ANYMORE: The length of fdm_sol must match the discretization, i.e. numel(opt_init.fdm_sol) = (opt_approx_method.n_int_1)*(opt_approx_method.n_int_2)*(system.dim)
        % NOW: If length of fdm_sol does not match discretization, fdm_sol is interpolated. Therefore: Carry out the additional checks in line 304 - 330
        % if isfield(opt_approx_method,'n_int_1') && isfield(opt_approx_method,'n_int_2')     % Use the provided value of n_int_1 and n_int_2 to check the size of fdm_sol      
        %     if length(opt_init.fdm_sol) ~= system.dim*opt_approx_method.n_int_1*opt_approx_method.n_int_2
        %         GC.error_msg{1,end+1} = append('The size of opt_init.fdm_sol is [', num2str(size(opt_init.fdm_sol,1)), ' x 1], but should be [(opt_approx_method.n_int_1)*(opt_approx_method.n_int_2)*(system.dim) x 1] = [', num2str(opt_approx_method.n_int_1*opt_approx_method.n_int_2*system.dim), ' x 1].');
        %         GC.error_msg{1,end+1} = 'You provided opt_approx_method.n_int_1 (n_int_1) and opt_approx_method.n_int_2 (n_int_2). If this error occurs, please identify n_int_1 and n_int_2 of your provided solution opt_init.fdm_sol and set n_int_1 as well as n_int_2 to their correct values.';
        %     end
        % elseif isfield(opt_approx_method,'n_int_1')                                         % Use the provided value of n_int_1 and the default value of n_int_2 to check the size of fdm_sol                   
        %     if length(opt_init.fdm_sol) ~= system.dim*opt_approx_method.n_int_1*50
        %         GC.error_msg{1,end+1} = append('The size of opt_init.fdm_sol is [', num2str(size(opt_init.fdm_sol,1)), ' x 1], but should be [(opt_approx_method.n_int_1)*(opt_approx_method.n_int_2)*(system.dim) x 1] = [', num2str(opt_approx_method.n_int_1*50*system.dim), ' x 1].');
        %         GC.error_msg{1,end+1} = 'You only provided opt_approx_method.n_int_1 (n_int_1). If this error occurs, please identify n_int_1 and n_int_2 of your provided solution opt_init.fdm_sol and set n_int_1 as well as n_int_2 to their correct values.';
        %     end
        % elseif isfield(opt_approx_method,'n_int_2')                                         % Use the default value of n_int_1 and the provided value of n_int_2 to check the size of fdm_sol                   
        %     if length(opt_init.fdm_sol) ~= system.dim*50*opt_approx_method.n_int_2
        %         GC.error_msg{1,end+1} = append('The size of opt_init.fdm_sol is [', num2str(size(opt_init.fdm_sol,1)), ' x 1], but should be [(opt_approx_method.n_int_1)*(opt_approx_method.n_int_2)*(system.dim) x 1] = [', num2str(50*opt_approx_method.n_int_2*system.dim), ' x 1].');
        %         GC.error_msg{1,end+1} = 'You only provided opt_approx_method.n_int_2 (n_int_2). If this error occurs, please identify n_int_1 and n_int_2 of your provided solution opt_init.fdm_sol and set n_int_1 as well as n_int_2 to their correct values.';
        %     end
        % else                                                                                % Use the default values of n_int_1 and n_int_2 to check the size of fdm_sol                   
        %     if length(opt_init.fdm_sol) ~= system.dim*50*50
        %         GC.error_msg{1,end+1} = append('The size of opt_init.fdm_sol is [', num2str(size(opt_init.fdm_sol,1)), ' x 1], but should be [(opt_approx_method.n_int_1)*(opt_approx_method.n_int_2)*(system.dim) x 1] = [', num2str(50*50*system.dim), ' x 1].');
        %         GC.error_msg{1,end+1} = 'You did not provide opt_approx_method.n_int_1 (n_int_1) and opt_approx_method.n_int_2 (n_int_2) and opt_init.fdm_sol does not match their default values.';
        %         GC.error_msg{1,end+1} = 'If this error occurs, please identify n_int_1 and n_int_2 of your provided solution opt_init.fdm_sol and set n_int_1 as well as n_int_2 to their correct values.';
        %     end
        % end
        GC.speak();
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

end