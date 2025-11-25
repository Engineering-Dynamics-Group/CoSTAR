% This local Gatekeeper methods checks the options structure for the solplot method.
%
% @obj:          Solution subclass object
% @DYN:          DynamicalSystem class object
% @options:      options structure for the solplot method
%
% THIS GATEKEEPER ONLY CHECKS THE FIELDs WHICH ARE NOT CHECKED IN SOLGET_GATEKEEPER 
% EXCEPTION: 'space' and 'resolution' needs to be checked, since this is needed inside this method

function options = solplot_gatekeeper(obj,DYN,options)

    GC = Gatekeeper();

    % Check if DYN is a Dynamical System class object
    GC.check_data(DYN,'DYN','DynamicalSystem',[],[]);    
    GC.speak('Error while using postprocessing method "solplot":');

    % Check if DYN and S have the same identification (if they are matching objects)
    if ~strcmpi(DYN.DYN_id,obj.S_id)
        GC.error_msg{1,end+1} = 'The DynamicalSystem class object and Solution class object have differing identification keys. Apparently you mixed objects and no postprocessing is done!';
    end
    
    % Check if option struct was generated using the costaropts function
    if ~isfield(options,'costaropts')
        GC.error_msg{1,end+1} = append('The contplot options structure was not created using the costaropts function, which is mandatory.');
    end

    % Check that solplot is not called for an equilibrium solution
    if strcmpi(DYN.sol_type,'equilibrium')
        GC.error_msg{1,end+1} = 'You are trying to use solplot for an equilibrium solution object, which is not meaningful.'; 
    end
    
    GC.speak('Error while using postprocessing method "solplot":');


    %% Check the fields
    options_mandatory_fieldnames  = {'zaxis','space'};      %mandatory fieldsnames in the options super structure. In 3D (e.g. hypertime quasi-periodic) you would always assume the third axis up, which is zaxis 
    options_allowed_fieldnames    = {'zaxis','space','xaxis','yaxis','index','mu','interval','resolution','figure','color','linestyle'};    %allowed fieldsnames in the options super structure

    allowed_space_values = {'time','trajectory','hypertime','frequency'};                                   %solution variables can in general not be plotted
    allowed_axis_values  = {'euclidean','all'};                                                             %allowed char values additional to function_handle for the axis
    allowed_colornames    = {'r','g','b','c','m','y','k'};                                                  %allowed char values for color definition
    allowed_linestyles    = {'-','--',':','-.'};                                                            %allowed char values for line style definition

    options     =  GC.fieldnames_to_lower(options);                                                         %Set all field names to lower case characters
    options     =  GC.fieldvalues_to_lower(options);                                                        %set fieldnames to lower case 

    GC.check_fields(options,'options',options_mandatory_fieldnames,options_allowed_fieldnames);            %updates the error_msg property of the gatekeeper        
    GC.speak('Error while using postprocessing method "solplot":');  


    %% Check the mandatory values -> 'zaxis' is checked for each solution space individually (see below)
    GC.check_data(options.space,'options.space','char', [] ,allowed_space_values);
    GC.speak('Error while using postprocessing method "solplot":');


    %% Check the optional field values  
    % -> 'interval', 'index' and 'mu' are checked in solget_gatekeeper
    % -> 'resolution' is checked here and in solget_gatekeeper. However, we need the check here for check_fcn_handle and the check needs to be BEFORE check_fcn_handle is called
    if isfield(options,'figure')
        GC.check_data(options.figure,'options.figure','matlab.ui.Figure',[],[]); 
        GC.speak('Error while using postprocessing method "solplot":');
    end 
   
    if isfield(options,'color')
        if isa(options.color,'char')
            GC.check_data(options.color,'options.color','char',[],allowed_colornames);
        elseif isa(options.color,'double')
            GC.check_data(options.color,'options.color','double',[1,3],[]);
            if max(options.color)>1||min(options.color)<0
                GC.error_msg{1,end+1} = append('You supplied an rgb vector via options.color, whose values must be in the range of [0,1]. However, your maximum is ',num2str(max(options.color)),' and your minimum is ',num2str(min(options.color)),'.');
            end
        else
            GC.error_msg{1,end+1} = append('The fieldvalue options.color of contplot options must either be a ''char'' or [1x3] array containing rgb values. Allowed ''char'' values are: ', GC.my_join(allowed_colornames));
        end
        GC.speak('Error while using postprocessing method "solplot":');
    end

    if isfield(options,'linestyle')
        GC.check_data(options.linestyle,'options.linestyle','char',[],allowed_linestyles);
        GC.speak('Error while using postprocessing method "solplot":');
    end

    if isfield(options,'resolution')
        % Check data type and dimension. Scalar AND [1x2] array only possible for quasi-periodic hypertime plots
        if strcmpi(DYN.sol_type,'quasiperiodic') && strcmpi(options.space,'hypertime')
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


    % Now check dependent on the solution space in which to plot
    switch options.space

        case 'time'

            options_time_mandatory_fieldnames  = {'zaxis','space'};                                                                        %mandatory fieldsnames in the options super structure
            options_time_allowed_fieldnames    = {'zaxis','space','index','mu','interval','resolution','figure','color','linestyle'};      %allowed fieldsnames in the options super structure
            
            GC.check_fields(options,'options',options_time_mandatory_fieldnames,options_time_allowed_fieldnames);                           %updates the error_msg property of the gatekeeper
            GC.speak('Error while using postprocessing method "solplot" with solution space ''time'':');
            
            GC.check_data(options.zaxis,'options.zaxis',{'function_handle','char'}, [] ,[]);
            if ~isa(options.zaxis,'function_handle') %zaxis is a string
                GC.check_data(options.zaxis,'options.zaxis','char',[],allowed_axis_values);
            else %zaxis is a function handle
                obj.check_fcn_handle(DYN,GC,options,'options.zaxis','matrix','matrix_or_vector');
            end
            GC.speak('Error while using postprocessing method "solplot" with solution space ''time'':');


        case 'trajectory'

            options_time_mandatory_fieldnames  = {'zaxis','space','xaxis'};                                                %mandatory fieldsnames in the options super structure
            options_time_allowed_fieldnames    = {'zaxis','space','xaxis','yaxis','index','mu','interval','resolution','figure','color','linestyle'};             %allowed fieldsnames in the options super structure

            GC.check_fields(options,'options',options_time_mandatory_fieldnames,options_time_allowed_fieldnames);            %updates the error_msg property of the gatekeeper
            GC.speak('Error while using postprocessing method "solplot" with solution space "trajectory":');
            
            GC.check_data(options.xaxis,'options.xaxis',{'function_handle'}, [] ,[]);
            GC.check_data(options.zaxis,'options.zaxis',{'function_handle'}, [] ,[]);
            if isfield(options,'yaxis');  GC.check_data(options.yaxis,'options.yaxis',{'function_handle'}, [] ,[]); end
            GC.speak('Error while using postprocessing method "solplot" with solution space ''trajectory'':');
            
            %All xxx_axis are now ensured to be function_handles
            obj.check_fcn_handle(DYN,GC,options,'options.zaxis','matrix','column_vector'); 
            obj.check_fcn_handle(DYN,GC,options,'options.xaxis','matrix','column_vector'); 
            if isfield(options,'yaxis');  obj.check_fcn_handle(DYN,GC,options,'options.yaxis','matrix','column_vector');  end     
            GC.speak('Error while using postprocessing method "solplot" with solution space ''trajectory'':');


        case 'hypertime'

            options_time_mandatory_fieldnames  = {'zaxis','space'};                                                 %mandatory fieldsnames in the options super structure
            options_time_allowed_fieldnames    = {'zaxis','space','figure','index','mu','resolution','color','linestyle'};      %allowed fieldsnames in the options super structure

            GC.check_fields(options,'options',options_time_mandatory_fieldnames,options_time_allowed_fieldnames);            %updates the error_msg property of the gatekeeper
            GC.speak('Error while using postprocessing method "solplot" with solution space ''hypertime'':');

            GC.check_data(options.zaxis,'options.zaxis',{'function_handle','char'}, [] ,[]);
            if ~strcmpi(class(options.zaxis),'function_handle')
                GC.check_data(options.zaxis,'options.zaxis','char',[],allowed_axis_values); 
            else %zaxis is a function handle
                obj.check_fcn_handle(DYN,GC,options,'options.zaxis','solution_argument','solution_argument'); 
            end
            
            if strcmpi(DYN.sol_type,'quasiperiodic') && isfield(options,'color')
                GC.error_msg{1,end+1} = 'The options field ''color'' is not allowed when plotting hypertime manifolds for quasi-periodic solutions.';
            end

            GC.speak('Error while using postprocessing method "solplot" with solution space ''hypertime'':');


        case 'frequency'
            
            options_time_mandatory_fieldnames  = {'zaxis','space'};                                                %mandatory fieldsnames in the options super structure
            options_time_allowed_fieldnames    = {'zaxis','space','index','mu','interval','resolution','figure','color','linestyle'};             %allowed fieldsnames in the options super structure
            
            GC.check_fields(options,'options',options_time_mandatory_fieldnames,options_time_allowed_fieldnames);            %updates the error_msg property of the gatekeeper
            GC.speak('Error while using postprocessing method "solplot" with solution space ''frequency'':');

            if isfield(options,'resolution')            % Make sure that 'resolution' is an even number (everything else has already been checked above)
                if mod(options.resolution,2) ~= 0       % ('resolution' is a scalar in solution space frequency)
                    GC.error_msg{1,end+1} = append('The data value of options.resolution is ', num2str(options.resolution), '. However, the resolution must be a positive even integer in solution space "frequency".');
                end
            end
            GC.speak('Error while using postprocessing method "solplot" with solution space ''frequency'':');

            GC.check_data(options.zaxis,'options.zaxis',{'function_handle','char'}, [] ,[]);
            if ~strcmpi(class(options.zaxis),'function_handle')
                GC.check_data(options.zaxis,'options.zaxis','char',[],allowed_axis_values);
            else %zaxis is a function handle
                obj.check_fcn_handle(DYN,GC,options,'options.zaxis','matrix','matrix_or_vector');
            end
            GC.speak('Error while using postprocessing method "solplot" with solution space ''frequency'':');

    end


    clear GC;


end