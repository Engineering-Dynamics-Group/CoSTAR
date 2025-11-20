%This local Gatekeeper methods checks the options structure for the contplot 
%method.
%
%@obj:          Solution subclass object
%@DYN:          DynamicalSystem class object
%@options:      options structure for the solplot method
%
%!!! THIS GATEKEEPER ONLY CHECKS THE FIELD WHICH ARE NOT CHECKED IN
%SOLGET_GATEKEEPER EXCEPTION: 'space' needs to be checked, since this is needed inside this method!!!

function options = contplot_gatekeeper(obj,DYN,options)

    GC = Gatekeeper();

    % Check if DYN is a Dynamical System class object
    GC.check_data(DYN,'DYN','DynamicalSystem',[],[]);
    GC.speak('You are using the SOLUTION.contplot method:');

    % Check if DYN and S have the same identification (if they are matching objects)
    if ~strcmpi(DYN.DYN_id,obj.S_id)
        GC.error_msg{1,end+1} = 'The DynamicalSystem class object and Solution class object have differing identification keys. Apparently you mixed objects and no postprocessing is done!';
    end

    % Check if option struct was generated using the costaropts function
    if ~isfield(options,'costaropts')
        GC.error_msg{1,end+1} = append('The contplot options structure was not created using the costaropts function, which is mandatory.');
    end

    GC.speak('Error while using postprocessing method "contplot":');
        

    %% Check the fields
    options_mandatory_fieldnames  = {'zaxis'};                                                      %mandatory fieldsnames in the options super structure
    options_allowed_fieldnames    = {'zaxis','figure','resolution','index','color','linestyle'};    %allowed fieldsnames in the options super structure
    
    allowed_axis_values  = {'min2','mean2','max2'};                                                         %allowed char values additional to function_handle for the axis
    allowed_colornames    = {'r','g','b','c','m','y','k'};                                                  %allowed char values for color definition
    allowed_linestyles    = {'-','--',':','-.'};                                                            %allowed char values for line style definition

    options     =  GC.fieldnames_to_lower(options);                                                         %Set all field names to lower case characters
    options     =  GC.fieldvalues_to_lower(options);                                                        %set fieldnames to lower case 

    GC.check_fields(options,'options',options_mandatory_fieldnames,options_allowed_fieldnames);            %updates the error_msg property of the gatekeeper        
    GC.speak('Error while using postprocessing method "contplot":');


           
    %% Check the optional field values  
    % -> 'resolution' is checked here and in solget_gatekeeper. However, we need the check here for check_fcn_handle and the check needs to be BEFORE check_fcn_handle is called
    % ->  'index' is checked in solget_gatekeeper, but there is an additional check for 'index' here
    if isfield(options,'figure')
        GC.check_data(options.figure,'options.figure','matlab.ui.Figure',[],[]); 
        GC.speak('Error while using postprocessing method "contplot":');
    end %else; options.figure = 1; 
    
    if isfield(options,'color')
        if isa(options.color,'char')
            GC.check_data(options.color,'options.color','char',[],allowed_colornames);
        elseif isa(options.color,'double')
            GC.check_data(options.color,'options.color','double',[1,3],[]);
            if max(options.color)>1||min(options.color)<0
                GC.error_msg{1,end+1} = append('You supplied an rgb vector via options.color, whose values must be in the range of [0,1]. However, your maximum  is ',num2str(max(options.color)),' and your minimum is ',num2str(min(options.color)),'.');
            end
        else
            GC.error_msg{1,end+1} = append('The fieldvalue options.color of contplot options must either be a ''char'' or [1x3] array containing rgb values. Allowed ''char'' values are: ', GC.my_join(allowed_colornames));
        end
        GC.speak('Error while using postprocessing method "contplot":');
    end

    if isfield(options,'linestyle')
        GC.check_data(options.linestyle,'options.linestyle','char',[],allowed_linestyles);
        GC.speak('Error while using postprocessing method "contplot":');
    end

    if isfield(options,'resolution')    % 'resolution' has to be checked here because it is needed for check_fcn_handle
        % Check data type and dimension. Scalar AND [1x2] array only possible for quasi-periodic hypertime plots
        if strcmpi(DYN.sol_type,'quasiperiodic')
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
        GC.speak('Error while using postprocessing method "contplot":');
    end 

    if isfield(options,'index')     % Additional check for 'index': Cannot be checked in solget_gatekeeper since this check is only required for contplot to work correctly 
        if ~ischar(options.index)
            if ~isempty(find(diff(options.index) ~= 1,1))   % Check that all elements are adjacent integers in increasing numbering, i.e. options.index(k+1) - options.index(k) = 1 for all k = 1,...,numel(options.index)-1
                GC.error_msg{1,end+1} = 'The elements of the options field "index" must be adjacent natural numbers in ascending order, e.g. the vectors [3,1,2] or [1:25,50:75] are not allowed.';
            end
        end
        GC.speak('Error while using postprocessing method "contplot":'); 
    end


    %% Check the mandatory values
    GC.check_data(options.zaxis,'options.zaxis',{'function_handle','char'}, [] ,[]);
    if ~isa(options.zaxis,'function_handle') %zaxis is a string
        GC.check_data(options.zaxis,'options.zaxis','char',[],allowed_axis_values);
    else %zaxis is a function handle
        obj.check_fcn_handle(DYN,GC,options,'options.zaxis','solution_argument','scalar');
    end

    GC.speak('Error while using postprocessing method "contplot":');


    clear GC;


end