%This function defines the input structs for costar in the right way.
%Bug in Matlab code:
%The way of defining a struct matters!
%
% Example:
%   csstruct = struct('a',{1,3},'b',2);
% leads to a different results than
%   csstruct.a =  {1,3};
%   csstruct.b =  2;
%
%This is especially true for using cells. This is problematic, since 
%the parameter key "param" is a cell.
%
% We NEED the second variant of defining a struct. This is ensured by 
% the costaropts function.
%
% The function call of costaropts is analogously to calling a struct
% 
% @varargin: Pairs of field names (char) and field values
% @csstruct: Output struct allocated in the right way (see above) with the additional field "costaropts".
%            This field is used to ensure, that the option structure was created using this function. It is
%            removed in the code later. 

function csstruct = costaropts(varargin)

    %% Check that there is an even number of input elements
    if size(varargin,1)>1
        small_gatekeeper('The input to costaropts must be a comma-separated list of field names and value pairs. ');
    end
    if mod(size(varargin,2),2)==1
        small_gatekeeper('Field names and value input arguments must come in pairs. There must be an even number of input arguments. ');
    end
    %% Check that there is no doubled entry
    tmp = cell(1,size(varargin,2)/2);
    counter = 0;
    for k = 1:2:numel(varargin)
        if isa(varargin{1,k},'char')
        counter = counter + 1;
        tmp{1,counter} = varargin{1,k};
         else
            small_gatekeeper(append('Your input argument',num2str(k),' should be of type "string" as it is a field name input but it is of type ', class(varargin{1,k})));
        end
    end
    
    %Check if there are any doubled entries
    [~,idx1,idx2] = unique(tmp);    
    if ~(numel(idx1)==numel(idx2))
        idx3 = setdiff(1:size(tmp,2),idx1);
        out = strjoin(tmp(1,idx3),' , ');
        small_gatekeeper(append('Your input options structure contains doubled entries for fieldname(s): ',out));
    end
    
 

    %Assigning the values to a struct
    for k = 1:2:numel(varargin)
            csstruct.(varargin{1,k}) = varargin{1,k+1};            
    end
    csstruct.costaropts = 'costaropts'; %This value is checked within costar, to make sure that the option structure was created using this function.

end


function small_gatekeeper(msg)
      disp('--------------------------------------------------');
      disp('---- The Gatekeeper says: You shall not pass! ----');
      disp('----  He has detected the following errors:  ----');
      disp(msg);
      disp('Example for correct input: costaropts("order",1,"rhs",Fcn,"dim",2).');
      error('Error in the costaropts function.');
end