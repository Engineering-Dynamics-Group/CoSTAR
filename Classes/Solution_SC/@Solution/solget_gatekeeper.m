% This function checks the options structure entries, which are common to solget, contget, solplot and contplont
%
% @obj:          Solution Class object
% @DYN:          DynamicalSystem Class object
% @options:      option structure for the postprocessing of the solution class

function options = solget_gatekeeper(obj,DYN,options)

    GC = Gatekeeper();

    %Check if DYN is a Dynamical System class object
    GC.check_data(DYN,'DYN','DynamicalSystem',[],[]);
    GC.speak('Error while using postprocessing method "solget":');

    %Check if DYN and S have the same identification (are matching objects)
    if ~strcmpi(DYN.DYN_id,obj.S_id)
        GC.error_msg{1,end+1} = 'The DynamicalSystem class object and Solution class object have differing identification keys. Apparently you mixed objects and no postprocessing is done!'; 
    end
    GC.speak('You are using the SOLUTION.solget/solplot/contplot method:'); 

    % Check if option struct was generated using the costaropts function
    if ~isfield(options,'costaropts')
        GC.error_msg{1,end+1} = append('The contplot options structure was not created using the costaropts function, which is mandatory.');
    end

    % Check that solget is not called for an equilibrium solution by the user (ATTENTION: solget can be called from contplot even for EQ, which is why options.call_from_contplot is needed)
    if isfield(options,'call_from_contplot')                    % If solget is called from contplot (which is fine)
        options = rmfield(options,'call_from_contplot');        % Remove the field (to not interfere with the allowed fieldnames) and proceed
    elseif strcmpi(DYN.sol_type,'equilibrium')                  % When solget is not called from contplot: Check if solution type is equilibrium. If yes: User called solget for EQ solution, which is not allowed
        GC.error_msg{1,end+1} = 'You are trying to use solget for an equilibrium solution object, which is not meaningful.';    % solget could also be called by solplot, but solplot is not allowed for EQ -> solget cannot be called by solplot for EQ
    end

    GC.speak('Error while using postprocessing method "solget":');


    %% Check the fields
    options_mandatory_fieldnames  = {'eval','space'};                                               %mandatory fieldsnames in the options super structure
    options_allowed_fieldnames    = {'eval','space','resolution','index','mu','interval'};          %allowed fieldsnames in the options super structure

    allowed_eval_values  = {'euclidean','all'}; %besides the function_handle call
    allowed_space_values = {'time','hypertime','frequency'};

    options = GC.fieldnames_to_lower(options);                                                         %Set all field names to lower case characters
    options = GC.fieldvalues_to_lower(options);                                                        %set fieldnames to lower case 

    GC.check_fields(options,'options',options_mandatory_fieldnames,options_allowed_fieldnames);            %updates the error_msg property of the gatekeeper        
    GC.speak('Error while using postprocessing method "solget":');


    %% Check the mandatory values
    
    GC.check_data(options.eval,'options.eval',{'function_handle','char'}, [] ,[]);
    if ~strcmpi(class(options.eval),'function_handle')
        GC.check_data(options.eval,'options.eval','char',[],allowed_eval_values); 
    end
    GC.check_data(options.space,'options.space','char', [] ,allowed_space_values);
    GC.speak('Error while using postprocessing method "solget":');


    %% Check the optional field values
    if isfield(options,'resolution')
        % Check data type and dimension. Scalar AND [1x2] array only possible for quasi-periodic hypertime plots using FDM currently
        if strcmpi(DYN.sol_type,'quasiperiodic') && strcmpi(options.space,'hypertime') && strcmpi(DYN.approx_method,'finite-difference')
            GC.check_data(options.resolution,'options.resolution','double',{'scalar','vector'},[]);
            if ~isscalar(options.resolution) && ~isequal(size(options.resolution),[1 2])            % If options.resolution is not a scalar: Check that is a [1x2] array
                GC.error_msg{1,end+1} = append('The data size of options.resolution is [', num2str(size(options.resolution)), ']. However, only scalars and [1x2] arrays are allowed.');
            end
        else
            GC.check_data(options.resolution,'options.resolution','double','scalar',[]); 
        end 
        % Check if value(s) are positive integer(s)
        if (any(mod(options.resolution,1)) || ~isempty(find(options.resolution<=0,1)))
            GC.error_msg{1,end+1} = append('The data value of options.resolution is [', num2str(options.resolution), ']. However, only positive integer value(s) are allowed.');
        end
        GC.speak('Error while using postprocessing method "solget", "solplot" or "contplot":');
    end 
    
    if isfield(options,'interval') 
        GC.check_data(options.interval,'options.interval','double',[1,2],[]); 
        if min(options.interval) < 0; GC.error_msg{1,end+1} = append('The minimum of options.interval is ',num2str(min(options.interval)),' but must be equal or larger to 0'); end
        if min(options.interval)>max(options.interval); GC.error_msg{1,end+1} = append('The minimum of options.interval is ',num2str(min(options.interval)),' is larger than its maximum ',num2str(min(options.interval)),', which is not allowed'); end    
        GC.speak('Error while using postprocessing method "solget" or "solplot":');
    end

    %Check that only mu or index has been set
    if isfield(options,'index') && isfield(options,'mu')
        GC.error_msg{1,end+1} = 'You provided options.index (the index of the continuation parameter) and options.mu (the continuation parameter itself). However, only either of the options is allowed.'; 
        GC.speak('Error while using postprocessing method "solget" or "solplot":'); 
    end

    if isfield(options,'index')     % If options.index is given 
        if ischar(options.index)
            GC.check_data(options.index,'options.index','char',[],'all');                   % options.index is allowed to have a numeric value or the string 'all'
        else
            GC.check_data(options.index,'options.index','double',{'scalar','vector'},[]);   % Instead of [], we could input "1:numel(obj.mu)". However, this leads to the Gatekeeper displaying all allowed values for the index (can be MANY numbers!)
            if any(mod(options.index,1))                                                    % Check that all elements are integers
                GC.error_msg{1,end+1} = 'There options field "index" contains non-integer number(s), which is not allowed.';
            end       
            if (min(options.index) < 1) || (max(options.index) > numel(obj.mu))             % Check that all elements of options.index are within the range of 1 to numel(obj.mu)
                GC.error_msg{1,end+1} = append('The value(s) of the options field "index" are in the range of [', num2str(min(options.index)), ',', num2str(max(options.index)), '], but they have to be in the range of [1,', num2str(max(numel(obj.mu))), '].');
            end
        end
        GC.speak('Error while using postprocessing method "solget", "solplot" or "contplot":'); 
    end

    if isfield(options,'mu')        % If options.mu is given 
        if ischar(options.mu)                                            
            GC.check_data(options.mu,'options.mu','char',[],'all');     %options.mu is allowed to have a numeric value or the string 'all'
            options.index = 1:size(obj.s,2);
        else
            GC.check_data(options.mu,'options.mu','double',{'scalar','vector'},[]); 
        end
        if (min(options.mu)<min(obj.mu))||(max(options.mu)>max(obj.mu))
            GC.error_msg{1,end+1} = ['min/max of options.mu is [',num2str(min(options.mu)),',',num2str(max(options.mu)),'] ','must be in the range of [',num2str(min(obj.mu)),',',num2str(max(obj.mu)),'].']; 
        end
        GC.speak('Error while using postprocessing method "solget" or "solplot":'); 
    end


    %Check for the individual spaces
    switch options.space 
        
        case 'time'
                time_mandatory_fieldnames  = {'eval','space'};                                                  %mandatory fieldsnames in the options super structure
                time_allowed_fieldnames    = {'eval','space','resolution','index','mu','interval'};             %allowed fieldsnames in the options super structure
                GC.check_fields(options,'options',time_mandatory_fieldnames,time_allowed_fieldnames);           %updates the error_msg property of the gatekeeper        
                GC.speak('Error while using postprocessing method "solget" with solution space "time":');

        case 'hypertime'
                hypertime_mandatory_fieldnames  = {'eval','space'};                                                 %mandatory fieldsnames in the options super structure
                hypertime_allowed_fieldnames    = {'eval','space','resolution','index','mu'};                       %allowed fieldsnames in the options super structure
                GC.check_fields(options,'options',hypertime_mandatory_fieldnames,hypertime_allowed_fieldnames);     %updates the error_msg property of the gatekeeper        
                GC.speak('Error while using postprocessing method "solget" with solution space "hypertime":');


        case 'frequency'
                frequency_mandatory_fieldnames  = {'eval','space'};                                                     %mandatory fieldsnames in the options super structure
                frequency_allowed_fieldnames    = {'eval','space','resolution','index','mu','interval'};                %allowed fieldsnames in the options super structure
                GC.check_fields(options,'options',frequency_mandatory_fieldnames,frequency_allowed_fieldnames);         %updates the error_msg property of the gatekeeper        
                GC.speak('Error while using postprocessing method "solget" with solution space "frequency":');
                if isfield(options,'resolution')            % Make sure that 'resolution' is an even number (everything else has already been checked above)
                    if mod(options.resolution,2) ~= 0       % ('resolution' is a scalar in solution space frequency)
                        GC.error_msg{1,end+1} = append('The data value of options.resolution is ', num2str(options.resolution), '. However, the resolution must be a positive even integer in solution space "frequency".');
                    end
                end
                GC.speak('Error while using postprocessing method "solget" or "solplot" with solution space "frequency":'); 

    end


    clear GC;   % Gatekeeper can get some rest for now


end