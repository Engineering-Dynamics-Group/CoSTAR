 % Check initial fields:
 %Check a option structure, whether all needed fields are there and not more than allowed. 
 %@data:                    struct or cell of strings to be checked
 %@name:                    string of the name of mystruct (needed for the error message)
 %@mandatory_fieldnames:    struct of strings: contains all fieldnames which must be contained in the struct to be checked
 %@allowed_fieldnames:      struct of strings: contains all fieldnames which are optional in the struct to be checked
 %@error_msg:               struct of error messages 

 %@varargin:                different variants of mandatory and allowed
                            %fields might be needed. These can be added by
                            %in varargin in: example:
                            %check_fields(obj,mystruct,name,mandatory1,allowed1,mandatory2,allowed2,mandatory3,allowed3).
                            %if mandatory1 and allowed1 fail, the next
                            %pairs are checked. If one combination of
                            %mandatory and allowed passes the test, no
                            %error is displayed

 function  check_fields(obj,data,name,mandatory_fieldnames,allowed_fieldnames,varargin)


        

         if isstruct(data)   %This is an option structure (as used by all sub-gatekeepers)
             %If check_fields is used within contplot_, solplot, solget_gatkeeper, this is necessary
             if isfield(data,'costaropts')
                 data = rmfield(data,'costaropts');
             end
             fields = fieldnames(data);
         else %This is a cell of strings (as used by the costarhelp class)
             fields = data;
         end
       

       mandatory{1,1}   = mandatory_fieldnames;
       allowed{1,1}     = allowed_fieldnames;
    %If there are multiple mandatory or allowed fields... fill them up. 

    if ~isempty(varargin)
        counter = 1;
        for k = 1:2:size(varargin,2)
            counter = counter + 1;
            mandatory{1,counter} = varargin{1,k};
            allowed{1,counter}   = varargin{1,k+1};
        end
    end


    %setdiff(A,B) returns the data in A that is not in B with no repetitions.
    %Is there all that I need? If setdiff is not empty, there are not all needed entries in B.
    %Is there more than allowed? If setdiff is not empty, there are more
    %fields than allowed.

    for k =1:size(mandatory,2)

        all_I_need{1,k} =        reshape(setdiff(mandatory{1,k},fields),1,[]);
        more_than_allowed{1,k} = reshape(setdiff(fields,allowed{1,k}),1,[]);
        
        if ~isempty(all_I_need{1,k});        passed(k,1) = 0; else passed(k,1) = 1; end
        if ~isempty(more_than_allowed{1,k}); passed(k,2) = 0; else passed(k,2) = 1; end
           
    end

    %Problem: both entries of mandatory and allowed need to pass the test
    %for a combination to be excepted. This is the job of passed.
    %This statement is only true, if there are somewhere two true entries
    %in a row.

    if ~any(passed(:,1).*passed(:,2))   %Is there a conflict in the fields?

            if size(mandatory,2)>1
                obj.error_msg{1,end+1} = append('The ', name,' structure has ', num2str(size(mandatory,2))',' sets of different field value options! There is a problem! Either:');
            end

        for k = 1:size(mandatory,2)     %Loop through all possibilities

            if ~isempty(all_I_need{1,k})    %Do I have all I need? If there is s.th. contained in the struct, this is what I need.
                temp{1,1} = append('"',all_I_need{1,k}{1,1},'"');
                for j =2:size(all_I_need{1,k},2)
                    temp{1,j} = append(', "',all_I_need{1,k}{1,j},'"');
                end

                if size(all_I_need{1,k},2)==1
                    obj.error_msg{1,end+1} = append('The ', name,' structure needs the field ',strjoin(temp),'.');
                else
                    obj.error_msg{1,end+1} = append('The ', name,' structure needs the fields ',strjoin(temp),'.');
                end
            end

            if ~isempty(more_than_allowed{1,k}) %Is there more than allowed? If there is s.th. contained in the struct, it is not allowed.

                temp{1,1} = append('"',more_than_allowed{1,k}{1,1},'"');
                for j =2:size(more_than_allowed{1,k},2)
                    temp{1,j} = append(', "',more_than_allowed{1,k}{1,j},'"');
                end

                if size(more_than_allowed{1,k},2)==1
                    obj.error_msg{1,end+1} = append('The ', name,' structure contains the field ',strjoin(temp),', which is not allowed. Allowed values are ', strjoin(allowed{1,k},' , '),'.');
                else
                    obj.error_msg{1,end+1} = append('The ', name,' structure contains the fields ',strjoin(temp),', which are not allowed. Allowed values are ', strjoin(allowed{1,k}),'.');
                end
            end

            if (size(mandatory,2)>1)&&(k<size(mandatory,2))
                obj.error_msg{1,end+1} = append('-------- OR --------');
            end



        end

    end

 end