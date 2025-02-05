% This function writes the log file during the computation
%
% Input arguments:
%   @DYN:       DynamicalSystem object
%   @varargin:  One or two char/string arguments. First argument of varargin is called @input and second argument is called @add_text
%       - @input:       General input text. There are special char/string values (see below) to control what to print
%       - @add_text:    Additional text, only possible for input = 'finalize' or input = 'finalize_error' 
%   
% Possible @input values:
%   - 'initialize':     Creates and sets up a text file when the computation starts
%   - 'diary_on':       Enables the diary function of MATLAB, which automatically prints command window output to the log file
%   - 'diary_off':      Disables the diary function
%   - 'finalize':       Prints the last information to the log file when the computation is finished
%   - 'finalize_error': Almost identical to 'finalize', but the footer message reports "Finished with error!"
%   - char/string:      General char/string values are directly printed into the log file

function write_log(DYN,varargin)

    %% Input handling

    % First argument "DYN" must be object of class DynamicalSystem
    if ~isa(DYN,'DynamicalSystem')
        error('Wrong data type of first input argument! Please make sure that it is an object of class DyncamicalSystem.')
    end

    % Two or three input arguments are possible
    switch nargin
        case 2
            input = varargin{1};
        case 3
            input = varargin{1};
            add_text = varargin{2};
            % Third argument "add_text" must be char or string
            if ~(isa(add_text,'char') || isa(add_text,'string'))   % Only char or string is allowed
                error('Wrong data type of second input argument! Please make sure that it is a char or a string.')
            end
            % Third argument "add_text" only allowed when input = 'finalize' or input = 'finalize_error'
            if ~(strcmpi(input,'finalize') || strcmpi(input,'finalize_error'))   % three
                error('Three input arguments are only allowed when ''finalize'' or ''finalize_error'' is specified as second input argument!')
            end
    end

    % Second argument "input" must be char or string
    if ~(isa(input,'char') || isa(input,'string'))
        error('Wrong data type of second input argument! Please make sure that it is a char or a string.')
    end


    %% Check whether a log file should be created
    if DYN.create_log == false
        return                                  % Leave function when no log file is desired
    end


    %% Definition
    filename = append('CoSTAR_Log_',DYN.DYN_id,'.txt');
    error_msg = 'Cannot open log file "%s". \nPlease do not touch the log file and do not change the current folder when CoSTAR is running.';


    %% Create or write the log file

    % Initialization
    if strcmpi(input,'initialize')

        fid = fopen(filename,'w');              % Open the log file with permission 'write' (new file is created)

        % Print header
        fprintf(fid,'%s\n%s\n%s\n%s\n%s\n',...
            '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', ...
            '%                                 CoSTAR                                  %', ...
            '%              Continuation of Solution Torus AppRoximations              %', ...
            append('%                              Version ',DYN.costar_version,'                                %'), ...
            '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'\n');
        fprintf(fid,'This is the log file of CoSTAR computation: %s\n',DYN.DYN_id);
        fprintf(fid,'\n');
        fprintf(fid,'---------------------------------------------------------------------------\n');
        fprintf(fid,'--------------------------------  Options  --------------------------------\n');
        fprintf(fid,'---------------------------------------------------------------------------\n');
        fprintf(fid,'\n');

        % Print the defined options
        option_structs = {'system','opt_sol','opt_approx_method','opt_init','opt_cont','opt_stability'};
        fprintf(fid,'CoSTAR was called using the following options:\n');
        fprintf(fid,'\n');
        for i = 1:numel(option_structs)
            if isempty(fieldnames(DYN.(option_structs{i})))                 % If user did not define an options struct
                str = '    (Not defined)';                                  % Set this text to print
                fprintf(fid,'options.%s:\n%s\n',option_structs{i},str);     % Print options struct
                fprintf(fid,'\n');                                          % Another new line is needed
            else    % If options struct was defined: Get the fields as formatted text (ATTENTION: needs version R2021a or above)
                str = formattedDisplayText(DYN.(option_structs{i}),'UseTrueFalseForLogical',true,'SuppressMarkup',true);
                fprintf(fid,'options.%s:\n%s\n',option_structs{i},str);     % Print options struct
            end
        end

        % Continue with initial solution
        fprintf(fid,'---------------------------------------------------------------------------\n');
        fprintf(fid,'---------------------------  Initial solution  ----------------------------\n');
        fprintf(fid,'---------------------------------------------------------------------------\n');
        fprintf(fid,'\n');
        fprintf(fid,'Starting computation at %s.\n',char(datetime('now','TimeZone','local','Format','yyyy-MMM-dd, HH:mm:ss')));
        fprintf(fid,'\n');
        tic;                                    % Start timer to record computation time
        
        fclose(fid);                            % Close the log file
    

    % Enable diary function
    elseif strcmpi(input,'diary_on')
        fid = fopen(filename,'a');              % Check if existing log file can be opened as required
        if fid == -1                            % Error handling
            error(error_msg,filename);
        end
        fclose(fid);                            % Close the log file (not needed for diary function)
        diary(filename)                         % Print recorded text to existing log file


    % Disable diary function
    elseif strcmpi(input,'diary_off')
        diary off


    % Finalization
    elseif strcmpi(input,'finalize') || strcmpi(input,'finalize_error')

        fid = fopen(filename,'a');              % Open the log file with permission 'append'
        if fid == -1                            % Error handling
            error(error_msg,filename);
        end

        % If there are three input arguments: Print additional text "add_text" first
        if nargin == 3                      
            fprintf(fid,'\n');
            fprintf(fid,append(add_text,'\n'));
        end

        % Print footer
        fprintf(fid,'\n');
        fprintf(fid,'---------------------------------------------------------------------------\n');
        if ~isempty(lastwarn)
            fprintf(fid,'-----------------------  Finished with warning(s)!  -----------------------\n');
        elseif strcmpi(input,'finalize_error')
            fprintf(fid,'-------------------------  Finished with error!  --------------------------\n');
        else
            fprintf(fid,'------------------------  Successfully finished!  -------------------------\n');
        end
        fprintf(fid,'---------------------------------------------------------------------------\n');
        fprintf(fid,'\n');
        fprintf(fid,'Computation finished at %s.\n',char(datetime('now','TimeZone','local','Format','yyyy-MMM-dd, HH:mm:ss')));
        fprintf(fid,'Elapsed time is %.3f seconds.',toc);

        fclose(fid);                            % Close the log file
        
        % Postprocessing: Remove printed links (links occur in fsolve output, but diary function prints them as ugly text using <link>)
        log_text = fileread(filename);                      % Read the whole log file
        log_text_new = eraseBetween(log_text,'<','>');      % Remove the links (everything between the '<' and '>')
        log_text_new = strrep(log_text_new,'<>','');        % Now remove the '<>', which are still there
        fid = fopen(filename,'w');                          % Open the log file with permission "write" and discard any existing content
        fwrite(fid,log_text_new);                           % Overwrite the log file with the adapted log file text
        fclose(fid);                                        % Close the log file


        
    % Write log file with arbitrary text
    else
        fid = fopen(filename,'a');              % Open the log file with permission 'append'
        if fid == -1                            % Error handling
            error(error_msg,filename);
        end
        fprintf(fid,append(input,'\n'));        % Print the text and jump to new line

        fclose(fid);                            % Close the log file


    end


end