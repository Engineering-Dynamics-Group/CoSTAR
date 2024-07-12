% This method is a static method from the costarhelp class
% It is called from the tutorial or minimal_example method
% Checks the input data and - if correct - searches for the corresponding files in the Tutorial folder
%
% @name:                String of the files we are looking for, either 'Tutorial' or 'Minimal_Example'
% @possible_words_1/2:  These are the possible word list for the user input (which are also used for calling the help functions)
% @alternative_1/2:     Alternative words for the possible_words_1/2! ATTENTION: They need to be equally long and in the same order.
%
% @found_file_names:    Contains the names of the files which correspond to the input data. The files can be opened with the edit function
% @other_file_names:    Contains the names of all other files which do not correspond to the input data, but still belong to @name


function   [found_file_names,other_file_names] = s_call_for_tutorial(name,possible_words_1,alternative_1,possible_words_2,alternative_2,varargin)

    varargin = varargin{1,1};

    key_words = {'',''};        % Initializing this cell with this voids makes the code for searching file containing the words below way easier
    all_words = {possible_words_1{1,:},alternative_1{1,:},possible_words_2{1,:},alternative_2{1,:}};

    
    %% Case: The user has specified the input
    if ~isempty(varargin)
        
        GC = Gatekeeper();
        
        % Check if the words contained in varargin are allowed
        if numel(varargin)>2; GC.error_msg{1,end+1} = append('There is a maximum of 2 allowed key words for the costarhelp. E.g. costarhelp.tutorial("key word 1","key word 2").'); end
        
            GC.check_fields(lower(varargin(1,:)),'key word for the costarhelp',{},lower(all_words)); %Are only allowed words in the user input?
            GC.speak('This is an error from the costarhelp function:');

        % Assign the words to key_words
        for k = 1:size(varargin,2)
            key_words{1,k} = varargin{1,k};
        end

        % Check if the two key_words are identical or from the same identifier group. If so: throw an error
        if size(varargin,2) == 2   %Only necessary for two words
            test = zeros(2,1);
            for k = 1:2 %We need only to check the first word group. If the word is not in the first word group - it is automatically in the second group (checked by Gatkeeper above)
                if any(strcmpi(key_words{1,k},{possible_words_1{1,:},alternative_1{1,:}}))
                    test(k,1) = 1;
                end
            end
            if prod(test)||prod(~test)
                GC.error_msg{1,end+1} = append('Your supplied key words ',key_words{1,1},' and ',key_words{1,2},' are either identical or from the same identifier group.');
                GC.error_msg{1,end+1} = append('Choose one key word from each identifier group.');
                GC.error_msg{1,end+1} = append('Identifier group 1 are: ',strjoin(possible_words_1,', '),' or alternatively ', strjoin(alternative_1,', '));
                GC.error_msg{1,end+1} = append('Identifier group 2 are: ',strjoin(possible_words_2,', '),' or alternatively ', strjoin(alternative_2,', '));
                GC.speak('This is an error from the costarhelp function:');
            end
        end

        clear GC;

        % Replace the key word with the expression from possible_words_x, if it is in the word group alternative_x
        for k = 1:size(varargin,2) %This is super ugly... please replace at some time
            
            tmp = strcmpi(key_words{1,k},alternative_1);                    %If it is in the alternative_1 group
            if any(tmp); key_words{1,k} = possible_words_1{1,tmp}; end      %Replace it by its counterpart
            
            tmp = strcmpi(key_words{1,k},possible_words_1);                 %If it is already in the lower caps possible_words_1 group... replace it by the correct spelling 
            if any(tmp); key_words{1,k} = possible_words_1{1,tmp}; end

            tmp = strcmpi(key_words{1,k},alternative_2);                    %see above
            if any(tmp); key_words{1,k} = possible_words_2{1,tmp}; end

            tmp = strcmpi(key_words{1,k},possible_words_2);                %see above
            if any(tmp); key_words{1,k} = possible_words_2{1,tmp}; end
      
        end

    end


    %% Find all tutorial files relevant for the requested case

    [current_path,~,~] = fileparts(mfilename('fullpath'));  % Get the current filepath

    Tutorials_path = append(current_path(1,1:end-19),'Tutorials');      % Get the Tutorials path. -19 removes 'Classes\@costarhelp' from current_path
    
    files = dir(fullfile(Tutorials_path, append(name,'_*','.m')));      % Get the names of all <name>_.m-files in the folder
                                                % files is a struct with fields name, folder, date, bytes, isdir, datenum
    n_found = 0;                                % Number of files found
    n_other = 0;                                % Number of other files

    files_found_tmp = cell(size(files));        % Initialize before loop
    other_files_tmp = cell(size(files));        % Initialize before loop

    for i = 1:size(files,1)                     % Loop for all files which have been found
        
        filename = split(files(i).name,'.m');   % Do not get "filename.m", but "filename" (needed for display output and function edit)
        if contains(filename{1,1},append('_',key_words{1,1})) && contains(filename{1,1},append('_',key_words{1,2}))     
            n_found = n_found + 1;
            files_found_tmp{n_found} = filename{1,1};       % All files that have been found according to the request are stored here
        else
            n_other = n_other + 1;
            other_files_tmp{n_other} = filename{1,1};       % If there are other files, they are stored here
        end 

    end

    found_file_names = files_found_tmp(1:n_found);
    other_file_names = other_files_tmp(1:n_other);


end