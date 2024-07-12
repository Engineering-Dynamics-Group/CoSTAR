% This function is a method of the costarhelp class and displays the 
% option help structure, defined in the corresponding s_help_opt_xxx files
% 
% help_struct:      Example: MUST contain the fields help_struct.n_fft.value and help_struct.n_fft.text

function s_disp_opt_help_struct(help_struct)


        names = fieldnames(help_struct);
        if ~isempty(names)
       
        max_width_prop = max([12,ceil(max(cellfun(@numel,names))*1.5)]); %12 is the length of "property"x1,5
        width_field_val  = 32;    
        width_string_val = 30;

        width_field_descri  = 48;    
        width_string_descri = 46;

        fprintf('%-*.*s%-*.*s%-*.*s\n', max_width_prop, max_width_prop, 'field name', width_field_val, width_string_val, 'allowed field values', width_field_descri, width_string_descri, 'description');   
        fprintf('%-s%-s%-s\n',    repmat('-',1,max_width_prop),repmat('-',1,width_field_val), repmat('-',1,width_field_descri));
        
        for k = 1:numel(names)
            sub_names = fieldnames(help_struct.(names{k,1})); %Get the field names of the struct within the struct
            
            %Get the formatted texts:
            value_text = costarhelp.s_format_string_to_cell(help_struct.(names{k,1}).(sub_names{1,1}),width_string_val,'character');        %Subdivide the single texts into a cell
            descri_text = costarhelp.s_format_string_to_cell(help_struct.(names{k,1}).(sub_names{2,1}),width_string_descri,'word');       %Subdivide the single texts into a cell

            cell_max_size = max(size(value_text,1),size(descri_text,1)); %Get the maximum number of cell entries size(,1) is necessary, since we need the first direction: value_text might only be a char and not a cell

            property_cell = vertcat(names{k,1},cell(cell_max_size-1,1));
            value_cell = vertcat(value_text,cell(cell_max_size-size(value_text,1),1));
            descri_cell = vertcat(descri_text,cell(cell_max_size-size(descri_text,1),1));
            
            for j = 1:cell_max_size
                fprintf('%-*.*s%-*.*s%-*.*s\n', max_width_prop, max_width_prop, property_cell{j,1}, width_field_val, width_string_val, value_cell{j,1}, width_field_descri, width_string_descri, descri_cell{j,1});
            end
            fprintf('\n');

        
        end
        end

end



