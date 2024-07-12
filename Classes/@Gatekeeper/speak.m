%Function for displaying error messages
%
%@obj:          Gatekeeper class object
%@varargin:     Additional cell of strings, which is displayed at the beginning of
%the message
%@error_msg: cell of string error messages

function speak(obj,varargin)
    
    if ~isempty(varargin)
        for k = 1:length(varargin)
            if ~ischar(varargin{1,k}); error('Error in Gatekeeper.speak: varargin passed to the function must be strings.'); end
        end
    end


    if ~isempty(obj.error_msg)
        disp('--------------------------------------------------');
        disp('---- The Gatekeeper says: You shall not pass! ----');
        disp('----   He has detected the following errors:  ----');

        for k = 1:length(varargin)
            disp(varargin{1,k});
        end
        for k=1:length(obj.error_msg)   %loop through the cell to display all messages
            disp(obj.error_msg{1,k});
        end
        disp('--------------------------------------------------');
        error('Gatekeeper error');
    end

end