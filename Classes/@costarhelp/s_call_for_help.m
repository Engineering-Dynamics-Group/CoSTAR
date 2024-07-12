%This method is a static method from the costarhelp class.
% It is called from the option structure help function. It 
% Checks the input data and - if correct - searches for the corresponding help files in the complete class folder
%
%
%@help_name:            String of the main help file we are looking for, e.g. opt_approx_method, opt_init
%@possible_words_1/2:   These are the possible word list for the user input (which are also used for calling the help functions)
%@alternative_1/2:      Alternative words for the possible_words_1/2! ATTENTION: They need to be equally long and in the same order.
%
%@help_method_names:    Contains the strings of the help_method_name(s), which correspond to the input data. The help_method_names string can be excecuted with evalc         
%@base_method_names:    Contains the strings of the base_method_name(s), which correspond to the help_name input. The base_method_name string can be excecuted with evalc         



function   [help_method_names,base_method_names] = s_call_for_help(help_name, possible_words_1,alternative_1,possible_words_2,alternative_2,varargin)

varargin = varargin{1,1};

key_words = {'',''};%Initializing this cell with this voids makes the code for searching file containing the words below way easier.
all_words = {possible_words_1{1,:},alternative_1{1,:},possible_words_2{1,:},alternative_2{1,:}};
%% Case: The user has specified the input
    if ~isempty(varargin)
        
        GC = Gatekeeper();
        %% Check if the words contained in varargin are allowed
        if numel(varargin)>2; GC.error_msg{1,end+1} = append('There is a maximum of 2 allowed key words for the costarhelp. E.g. costarhelp.opt_approx_method("key word 1","key word 2").'); end
        
            GC.check_fields(lower(varargin(1,:)),'key word for the costarhelp',{},lower(all_words)); %Are only allowed words in the user input?
            GC.speak('This is an error from the costarhelp function:');

        % Assign the words to key_words
        for k = 1:size(varargin,2)
            key_words{1,k} = varargin{1,k};
        end

        %% Check if the two key_words are identical or from the same identifier group. If so: throw an error
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
        %% Replace the key word with the expression from possible_words_x, if it is in the word group alternative_x

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



    %Find all help files relevant for the requested case:

     %% Select all help files: they begin all with s_help*.m 

    [tmp_file,~,~] = fileparts(mfilename('fullpath'));   %Get the current path of the @costarhelp folder
    class_filepath = tmp_file(1,1:end-12); %This is the filepath of the class folder: "-12" subtracts the "\@costarhelp" part of the filepath

    d = dir([class_filepath,'/**/*.*']); %Get all folders and subfolders of the classfolder (that is what '**/*.*' is for)
    dfolders = d([d(:).isdir]); % remove all files (isdir property is 0)
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'})); % remove '.' and '..' 

    %Create a default table
    help_files = array2table(zeros(0,2)); 
    help_files.Properties.VariableNames = {'file','folder'};
    
    base_files = array2table(zeros(0,2)); 
    base_files.Properties.VariableNames = {'file','folder'};    %These are all files, which contain only the option string


    %Loop through all folders and subfolders: Filter all files starting with s_help_help_name and then only keep the files with none, one or two keywords. 
    for k = 1:numel(dfolders)
        folder = fullfile(dfolders(k,1).folder,dfolders(k,1).name);
        tmp_file = dir(fullfile(folder,append('s_help_',help_name,'_*.m'))); %Is there a m-file starting with s_help_"help_name"?
        
        if ~isempty(tmp_file)
            
            if contains(tmp_file.name,append('_',key_words{1,1}))&&contains(tmp_file.name,append('_',key_words{1,2})) %This works always: If no key word was supplied by the user, we are simply searching for scripts with an underscore "_"
                help_files = [help_files;{tmp_file.name,tmp_file.folder}];  %Then add the folder name and the folder path to the default table
            else
                base_files = [base_files;{tmp_file.name,tmp_file.folder}];
            end
        end 
    end

    %% Build up executable file names 
   help_method_names = cell(size(help_files,1),1);
   for k = 1:size(help_files,1)

        method_name = split(help_files(k,:).file,'.m');   %Do not get "test_1.m", but "test_1" (needed for evalc to work)
        class_name =  split(help_files(k,:).folder,'@');


        help_method_names{k,1} = append(class_name{end,1},'.',method_name{1,1},'()');
        help_method_names{k,2} = class_name{end,1};
        help_method_names{k,3} = method_name{1,1};

   end


   base_method_names = cell(size(base_files,1),1);
   for k = 1:size(base_files,1)

        method_name = split(base_files(k,:).file,'.m');   %Do not get "test_1.m", but "test_1" (needed for evalc to work)
        class_name =  split(base_files(k,:).folder,'@');


        base_method_names{k,1} = append(class_name{end,1},'.',method_name{1,1},'()');
        base_method_names{k,2} = class_name{end,1};
        base_method_names{k,3} = method_name{1,1};

   end
    

end




% costarhelp.opt_approx_method('FGM','PS');