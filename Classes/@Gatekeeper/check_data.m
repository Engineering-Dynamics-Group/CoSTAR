%Checks option structure entries for the correct type, the allowed values and the allowed value size.
%However, you do not need to supply all values. If only [] is supplied...
%the check is not done.
%If an error occurs, the error_msg property of the Gatekeeper is updated.
%
%@obj:      Gatekeeper object
%@name:     string for potential error message output
%@data:     data to be checked 
%@data_type:     allowed data type as a char or a cell of chars, e.g.
%'double' or {'double','function_handle'}
%@data_value:    allowed data values, also the strings "positive" or "negative" are allowed to indicate that the value must strict  
%@data_size:     allowed data size values are an array e.g. [1,10] or the
%strings 'scalar','vector','matrix' or a cell of a combination of these
%e.g. {'scalar','vector'}


function check_data(obj,data,name,data_type,data_size,data_value)


    %Check the usage of the function itself:
    if ~isempty(data_size)
        if isa(data_size,'double')
            key = 'explicit';
            if (~(size(data_size,1)==1))&&(~(size(data_size,2)==2)); error('Gatekeeper - check_data: data_size must be a 1x2 array.'); end
        elseif isa(data_size,'char')||isa(data_size,'cell')
            key = 'general';
            if (~(size(data_size,1)==1))&&(~isempty(setdiff(data_size,{'scalar','vector','matrix'}))); error('Gatekeeper - check_data: data_size must be scalar, vector or matrix.'); end
        else
            error('Gatekeeper - check_data: data_size must either be a 1x2 array or the chars scalar, vector or matrix.');
        end
    end

    %% Now do what the function is supposed to do

    if isempty(data)    
            obj.error_msg{1,end+1} = append('The data field ',name,' is empty. Please supply a value.'); 
    else
        %Check data type
        if (~isempty(data_type))&&(~isempty(setdiff({class(data)},data_type)))
            obj.error_msg{1,end+1} = append('The data type of ',name,' is of type "',class(data),'" but should be "',obj.my_join(data_type),'".'); 
        end
        obj.speak();

        %Check data size (upper statement made sure, that there is no other possibility for data_size to be an array or a char)
        if (~isempty(data_size))
            if strcmpi(key,'explicit')      
                if (~all(size(data)==data_size))
                    obj.error_msg{1,end+1} = append('The data size of ',name,' is ',obj.my_join(size(data)),' but should be ',obj.my_join(data_size),'.'); 
                end
            else 
                if numel(size(data))==2 %data is a [n x m] array
                    tmp = 'none';
                    if length(data) == 1; tmp = 'scalar'; elseif (length(data)>1)&&(length(data)==numel(data)); tmp = 'vector'; else tmp = 'matrix'; end %Determine data size: scalar, vector, matrix
                    if ~isempty(setdiff(tmp,data_size)); obj.error_msg{1,end+1} = append('The data size type of ',name,' is ',tmp,' but should be ',obj.my_join(data_size),'.'); end %IF data size type is not equal to requested... throw error         
                else
                    obj.error_msg{1,end+1} = append('The data size of ',name,' must be a [n x m] array.'); 
                end
            end
        end
        obj.speak();        
        
        %Data Value is never allowed to contain the values NaN or Inf (which is only sensible, if data is numeric) or be complex
        if isnumeric(data)
            if ~all(~isnan(data),'all')
                obj.error_msg{1,end+1} = append('The data value of ',name,' contains "NaN" values, which is not allowed.');
            end
            if ~all(~isinf(data),'all')
                obj.error_msg{1,end+1} = append('The data value of ',name,' contains "Inf" values, which is not allowed.');
            end
            if ~isreal(data)
                obj.error_msg{1,end+1} = append('The data value of ',name,' is ',obj.my_join(data),' and contains imaginary values, which is not allowed.');
            end
        end
        obj.speak();

        %Data_value accepts the strings 'positive' and 'negative', which indicate that all the data must be strictly positive or negative.
        if ~isempty(data_value)&&(ischar(data_value)||iscellstr(data_value))
            if ischar(data_value)
                if (strcmpi(data_value,'positive')||strcmpi(data_value,'negative'))
                    if strcmpi(data_value,'positive')
                        if any(data<0)
                            obj.error_msg{1,end+1} = append('The data value of ',name,' is ',obj.my_join(data),'. However, only positive values are allowed.');
                        end
                    elseif strcmpi(data_value,'negative')
                        if any(data>0)
                            obj.error_msg{1,end+1} = append('The data value of ',name,' is ',obj.my_join(data),'. However, only negative values are allowed.');
                        end

                    end
                end
            elseif(~isempty(setdiff(data,data_value)))
                obj.error_msg{1,end+1} = append('The data value of ',name,' is ',obj.my_join(data),'. Allowed values are: ',obj.my_join(data_value),'.'); 
            end
        elseif (~isempty(data_value))&&isnumeric(data_value)
               if (~isempty(setdiff(data,data_value)))
                obj.error_msg{1,end+1} = append('The data value of ',name,' is ',obj.my_join(data),'. Allowed values are: ',obj.my_join(data_value),'.'); 
               end
        end
        obj.speak();
      

    end

end