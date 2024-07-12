% This function calls all files in the Tutorial folder
% All figures, displayed messages, warnings or errors are suppressed.
% Files are allowed to contain any "clear all" statements 
%
% @varargin{1,1}:   Give the solution type EQ, PS or QPS
% @varargin{1,2}:   Give the approximation method FDM, FGM or SHM

function check_tutorials(varargin)
  
    % addpath(genpath('..\'))                               % Add main folder of CoSTAR and all subfolders to search path
    clearvars -except varargin



    %% Get the file names of the Tutorials folder
  
    [current_path,~,~] = fileparts(mfilename('fullpath'));  % Get the current filepath

    cd(current_path)
        cd ..\Tutorials
        Tutorials_path = cd;    %                           % Get the Tutorials path
    cd(current_path)
    
    files = dir(fullfile(Tutorials_path, '*.m'));           % Get the names of all .m-files in the folder
   
  

    %% Get the files to be executed

    code_word = {'',''};

    if nargin == 0                    
        % Nothing to do here -> all files are to be checked
    elseif nargin == 1
        code_word{1,1} = varargin{1,1};
    elseif nargin == 2
        code_word{1,1} = varargin{1,1};
        code_word{1,2} = varargin{1,2};
    end
    

    n_files = 0;                                    % Number of files to be executed
    filenames = cell(size(files,1),1);              % Initialize before loop

    for i = 1:size(files,1)                         % Loop for all files which have been found
        
        filename = split(files(i).name,'.m');       % Do not get "filename.m", but "filename" (needed for evalc to work)
        
        % If filename contains '_<code_word_1>' and '_<code_word_2>'
        if contains(filename{1,1},append('_',code_word{1,1})) && contains(filename{1,1},append('_',code_word{1,2}))
            n_files = n_files + 1;
            filenames{n_files} = filename{1,1};     % filename is a [2x1] array, but filename{2,1} is empty
        end 

    end

    save('files.mat','n_files','filenames')         % Save important variables because "clear" in the tutorial files clears all variables



    %% Do the testing
    set(0,'DefaultFigureVisible','off');            % Suppress all figure windows
    set(0,'DefaultFigureWindowStyle','docked');

    message = cell(n_files,2);                      % Array storing the output messages

    for i = 1:n_files

        save('temp.mat','i','message')              % Save important variables because "clear" in the tutorial files clears all variables

        try

            lastwarn('');                               % Prepare for checking warning. Clear last warning message
            output = evalc(filenames{i});               % evalc catches any output by the function and does not display it
            clearvars; close all;                       % Close all figures if some figures have been opened
                                               
            load('temp.mat','i','message')              % Load missing variables
            load('files.mat','n_files','filenames')     % Load missing variables

            [warnMsg, ~] = lastwarn;                    % Get the last warning message
            if ~isempty(warnMsg)
                message{i,1} = append('warning: ',warnMsg);
                message{i,2} = 'warning';
            else
                message{i,1} = 'passed';
                message{i,2} = 'passed';
            end

        catch ME

            load('temp.mat','i','message')              % Load missing variables
            load('files.mat','n_files','filenames')     % Load missing variables
            message{i,1} = append('error: ',ME.message);
            message{i,2} = 'error';

        end

    end
    

    % Create the display text (has been deleted as well!)
    for i = 1:n_files
        cprintf('text','%s',filenames{i});
        if strcmpi(message{i,2},'passed')
            cprintf('green',' %s \n','passed');                 % green
        elseif strcmpi(message{i,2},'warning')
            cprintf('SystemCommands',' %s \n',message{i,1});    % orange
        elseif strcmpi(message{i,2},'error')
            cprintf('Errors',' %s \n',message{i,1});            % red
        end
    end


    delete('files.mat')                                 % Delete file - not required anymore
    delete('temp.mat')                                  % Delete file - not required anymore


    set(0,'DefaultFigureVisible','on');                 % Allow all figure windows again
    set(0,'DefaultFigureWindowStyle','normal')


    beep;                                               % Produce operating system beep sound
    


end