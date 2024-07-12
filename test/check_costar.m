%This function calls all test files in this folder, which start with
%"test_xxx" (underscore is relevant). This way you can automatically
%test all kinds of different test function from different developers.
%
%CoSTAR_test suppresses all figures, displayed messages, warnings or
%errors.
%
% Test files are now allowed to contain any "clear all;" statements 
%
%@varargin{1,1}:    give the solution type EQ,PS or QPS or the approximation method Shooting, FGM, etc. which is then exclusively evaluated
% @varargin{1,2}:   give the function the path to the costar version to
%be checked

function check_costar(varargin)
  
    clearvars -except varargin

    %% Add the path to the costar version of the current folder or the path to the costar version given in varargin to the system path
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(genpath(filepath));                         %So that the data folder is also included.
   
    select = 0;
    cd(filepath)
        cd .. 
        code_path = cd;
    cd(filepath)

    code_word_2 = [];
    select_1 = 0;
    select_2 = 0;

    addpath(genpath(code_path));
    if isempty(varargin)
        %nothing to do here
    elseif size(varargin,2) == 1
        code_word_1 = varargin{1,1};
        select_1 = 1;
    elseif size(varargin,2) == 2
        code_word_1 = varargin{1,1};
        code_word_2 = varargin{1,2};
        select_2 = 1;
    end
  

    %% Select the files
    files = dir(fullfile(filepath, 'test_*.m')); %Get the names of all .m-files in the folder, which start with test_
    counter = 0;
    for ii = 1:size(files,1)

        if select_1
            tmp0 = split(files(ii).name,'.m');   %Do not get "test_1.m", but "test_1" (needed for evalc to work)
            if contains(tmp0{1,1},append('_',code_word_1))
                counter = counter + 1;
                tmp{1,counter} = tmp0{1,1};
            end   
        elseif select_2
            tmp0 = split(files(ii).name,'.m');   %Do not get "test_1.m", but "test_1" (needed for evalc to work)
            if contains(tmp0{1,1},append('_',code_word_1))&&contains(tmp0{1,1},append('_',code_word_2))
                counter = counter + 1;
                tmp{1,counter} = tmp0{1,1};
            end 

        else
            tmp0 = split(files(ii).name,'.m');   %Do not get "test_1.m", but "test_1" (needed for evalc to work)
            tmp{1,ii} = tmp0{1,1};
        end

    end

    %% Do the testing
    set(0,'DefaultFigureVisible','off'); %Suppress all figure windows
    set(0,'DefaultFigureWindowStyle','docked');

    for ii = 1:size(tmp,2)

        script = tmp{1,ii};
        cprintf('text','%s',script);
    
        try
            lastwarn(''); %Prepare for checking warning. Clear last warning message
            tmp0 = evalc(script);    %evalc catches any output by the function and does not display it
            clearvars -except ii files tmp script
            close all; %If some figure get opened anyway
    
            [warnMsg, ~] = lastwarn;
            if ~isempty(warnMsg)
                cprintf('SystemCommands',' %s \n',append('warning: ',warnMsg));   %orange
            else
                cprintf('green',' %s \n','passed'); %green
            end
    
        catch ME
            cprintf('Errors',' %s \n',append('error: ',ME.message));    %red
        end
    
    
    
    end
    
    set(0,'DefaultFigureVisible','on'); %Allow all figure windows
    set(0,'DefaultFigureWindowStyle','normal')

    beep;
    


end