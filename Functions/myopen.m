
% Help function for the development crew... it opens all m files containing the "string" at the beginning

function myopen(string)


    %% Select all help files: they begin all with s_help*.m 
 
    [tmp_file,~,~] = fileparts(mfilename('fullpath'));   %Get the current path of the @costarhelp folder
    costar_filepath = tmp_file(1,1:end-9); %This is the filepath of the class folder: "-9" subtracts the "function" part of the filepath

    d = dir([costar_filepath,'/**/*.*']); %Get all folders and subfolders of the classfolder (that is what '**/*.*' is for)
    dfolders = d([d(:).isdir]); % remove all files (isdir property is 0)
    dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'})); % remove '.' and '..' 

        for k = 1:numel(dfolders)
            folder = fullfile(dfolders(k,1).folder,dfolders(k,1).name);
            tmp_file = dir(fullfile(folder,append('*',string,'*.m'))); %Is there a m-file starting with s_help_"help_name"?

            if ~isempty(tmp_file)

                for i = 1:length(tmp_file)
                    if contains(tmp_file(i).name,string) %This statement is Case sensitive in contrast to dir
                        open(fullfile(tmp_file(i).folder,tmp_file(i).name));
                    end
                end
             
            end

        end



end