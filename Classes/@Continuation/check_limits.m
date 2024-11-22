% Function checks whether calculated curve point is in the boundaries given
% by mu_limit and the maximal number of iterations is not surpassed. It
% furthermore displays information in the command window
%
% @obj:      Continuation class object

function obj = check_limits(obj,DYN)

obj.p_local_cont_counter = obj.p_local_cont_counter + 1;                % Set counter to next integer
local_cont_counter = obj.p_local_cont_counter;

formatSpec = '%.4f';

% Create info text and log file entry
info_text = append('Iter: ',num2str(local_cont_counter),' -- mu = ',num2str(obj.p_mu0,formatSpec),...
    ' -- stepwidth = ',num2str(obj.step_width,formatSpec));
if strcmpi(DYN.stability,'on')
    if obj.p_n_unstable_1 == 0;     info_text = append(info_text,' -- stable');
    elseif obj.p_n_unstable_1 > 0;  info_text = append(info_text,' -- unstable');   end
end
write_log(DYN,info_text)

% Display continuation step information
if strcmpi(DYN.display,'off')
    % Display nothing
elseif strcmpi(DYN.display,'final') || strcmpi(DYN.display,'iter')
    % The info text is only displayed for the latest N solutions in the command window
    % Idea: Save the info text of the latest N solutions in obj.p_last_msg. 
    % In the next loop, delete these info texts from the command window, update obj.p_last_msg and display it
    % Furthermore, warning messages are not deleted, but instead iteratively passed to the top (they are printed again at the top to have the orange color)
    %  (i.e. the warning messages are gathered at the top and N info texts are displayed below)
    if strcmpi(DYN.display,'final');    N = 1;                          % Display only the latest info text in case of DYN.display = 'final'
    elseif strcmpi(DYN.display,'iter'); N = 5; end                      % Display the latest 5 info texts in case of DYN.display = 'iter'
    reverse_str = repmat(sprintf('\b'), 1, length(obj.p_last_msg));     % This string will delete the old displayed info texts from the previous loop
    obj.p_last_msg = sprintf('%s%s\n',obj.p_last_msg,info_text);        % Add the new info text
    if local_cont_counter > (N+1)                                       % Remove first line when at least N info texts have been displayed (+1 for initial solution)
        end_str = append('Iter: ',num2str(local_cont_counter-N+1));     % String identifying where to end the deletion
        delete_str = extractBefore(obj.p_last_msg,end_str);             % This is the part which needs to be deleted
        if startsWith(obj.p_last_msg,'Warning')                         % If delete_str begins with a warning message (i.e. the warning message is the top line)
            delete_str_split = splitlines(delete_str(1:end-1));         % Split delete_str in its lines, because there can be multiple warnings
        elseif contains(delete_str,'Warning')                           % If delete_str contains a warning message, but there is an info text first
            delete_str = extractBefore(delete_str,'Warning');           % Exclude the warning message from deletion and delete only the info text (i.e. the top line)
        end
        obj.p_last_msg = erase(obj.p_last_msg,delete_str);              % Remove the top info text
    end
    if exist('delete_str_split','var')                                  % If we extracted a warning (i.e. a warning has been passed to the top)
        fprintf(reverse_str);                                           % Remove the info texts and the warning,
        for i = 1:numel(delete_str_split)                               % Scan the splitted delete_str for warnings
            if strcmpi(delete_str_split{i}(1:8),'Warning:'); warning(delete_str_split{i}(10:end)); end      % Print the warning(s) again to have the orange color
        end  
        fprintf([obj.p_last_msg]);                                      % Now print the new info texts
    else                                                                % Normal case: No warning has been passed to the top
        fprintf([reverse_str,obj.p_last_msg]);                          % Display the latest N info texts in command window (including warnings)
    end                             
else 
    % In case of DYN.display \in {'iter-detailed','step-control','error-control','full'}, the info text is displayed for each solution
    disp(info_text)
end


% Check the stopping criteria
if(obj.p_mu0>obj.mu_limit(1,2))         % Check if calculated curve point is above upper limit of mu
    obj.p_contDo = 0;
    clear fsolve;
    stopping_text = append('CoSTAR stopped because maximal value of continuation parameter mu_max = ', num2str(obj.mu_limit(1,2)), ' was reached.');
    write_log(DYN,'finalize',stopping_text)
    if ~strcmpi(DYN.display,'off');  disp(' '); disp(stopping_text);  end
    obj.p_stopping_flag = stopping_text;

elseif(obj.p_mu0<obj.mu_limit(1,1))     % Check if calculated curve point is below bottom limit of mu
    obj.p_contDo = 0;
    clear fsolve;
    stopping_text = append('CoSTAR stopped because minimal value of continuation parameter mu_min = ', num2str(obj.mu_limit(1,1)), ' was reached.');
    write_log(DYN,'finalize',stopping_text)
    if ~strcmpi(DYN.display,'off');  disp(' '); disp(stopping_text);  end
    obj.p_stopping_flag = stopping_text;

elseif(obj.p_local_cont_counter>=obj.max_cont_step)      % Check if maximal number of continuation points has been surpassed
    obj.p_contDo = 0;
    clear fsolve;
    stopping_text = append('CoSTAR stopped because maximal number of continuation steps max_cont_step = ', num2str(obj.max_cont_step), ' was reached.');
    write_log(DYN,'finalize',stopping_text)
    if ~strcmpi(DYN.display,'off');  disp(' '); disp(stopping_text);  end
    obj.p_stopping_flag = stopping_text;

else
    write_log(DYN,'\n-------------------------------------')

end

end