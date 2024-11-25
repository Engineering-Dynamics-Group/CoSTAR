% This method is a a submethod for the SOLUTION class gatekeepers to check the function handles
%
% @obj:     SOLUTION subclass object
% @DYN:     DYNAMICALSYSTEM class object
% @GC:      GATEKEEPER class object
% @fcn_handle_options:  can be a function_handle to be checked OR an options struct with field 'eval' or 'zaxis' containing a function handle
% @name_handle:         name of the function_handle as 'char' for error message
% @d_in:    dimension of input of the fcn_handle: can be
%            - char: 'solution_argument', 'matrix', 'matrix_or_vector', 'column_vector', 'row_vector' or 'scalar'.
%                     -> 'solution_argument' automatically selects the right dimension for an argument of the corresponding solution type
%            - an array specifying the dimension
% @d_out:   excpected dimension of the output of the function handle: can be
%            - char: 'solution_argument', 'matrix', 'matrix_or_vector', 'column_vector', 'row_vector' or 'scalar'.
%                     -> 'solution_argument' automatically selects the right dimension for an argument of the corresponding solution type
%            - an array specifying the dimension

function check_fcn_handle(obj,DYN,GC,fcn_handle_options,name_handle,d_in,d_out)


%% Get parameters

% Get the function handle
if isa(fcn_handle_options,'function_handle')    
    fcn_handle = fcn_handle_options;            % fcn_handle_options is a function handle
elseif isfield(fcn_handle_options,'eval')       
    fcn_handle = fcn_handle_options.eval;       % This is the case for the call from solget_gatekeeper
elseif isfield(fcn_handle_options,'zaxis')      
    fcn_handle = fcn_handle_options.zaxis;      % This is the case for the call from contplot_gatekeeper and solplot_gatekeeper
end


% Get the test sizes
dim = DYN.dim;                                  % To test fcn_handle output array dimension associcated with the state variables

if isfield(fcn_handle_options,'resolution')
    test_res = fcn_handle_options.resolution;   % To test fcn_handle output array dimension(s) associcated with (hyper)time direction
elseif strcmpi(DYN.sol_type,'quasiperiodic') && ~isfield(fcn_handle_options,'space')    % This is for the case when contplot is called for a quasi-periodic solution
    test_res = 50;                              % 50 is default size for quasi-periodic hypertime space
elseif strcmpi(DYN.sol_type,'quasiperiodic') && isfield(fcn_handle_options,'space') && strcmpi(fcn_handle_options.space,'hypertime')    % This is for the case when solplot is called for a quasi-periodic solution
    test_res = 50;                              % 50 is default size for quasi-periodic hypertime space
else
    test_res = 200;                             % 200 is the default size for everything else
end
if isfield(fcn_handle_options,'space') && strcmpi(fcn_handle_options.space,'frequency')
    test_res = test_res/2;                      % We need to half the test size for frequency solution space due to the FFT
end



%% Set the input and output array sizes

% Input array size
if isa(d_in,'char')

    switch d_in
        case 'solution_argument'
            switch DYN.sol_type
                case 'equilibrium'
                    d_in = [dim,1];
                case 'periodic'
                    d_in = [test_res,dim];
                case 'quasiperiodic'
                    d_in = [test_res,test_res,dim];
                otherwise
                    error('Something went wrong in SOLUTION class method check_fcn_handle.');
            end

        case 'matrix'
            d_in = [test_res,dim];

        case 'column_vector'
            d_in = [test_res,1];

        case 'row_vector'
            d_in = [1,test_res];

        case 'scalar'
            d_in = [1,1];

        otherwise
            error('Something went wrong in SOLUTION class method check_fcn_handle.');

    end

end


% Output array size
if isa(d_out,'char')

    switch d_out
        case 'solution_argument'
            switch DYN.sol_type
                case 'equilibrium'
                    d_out = [1,1];
                case 'periodic'
                    d_out = [test_res.*ones(dim,1), (1:dim)'];
                case 'quasiperiodic'
                    d_out = [test_res.*ones(dim,1), test_res.*ones(dim,1), (1:dim)'];
                otherwise
                    error('Something went wrong in SOLUTION class method check_fcn_handle.');
            end

        case 'matrix'
            d_out = [test_res,test_res];  % test_size is a random number

        case 'matrix_or_vector'
            d_out = [test_res.*ones(dim,1), (1:dim)'];

        case 'column_vector'
            d_out = [test_res,1];          % test_size is a random number

        case 'row_vector'
            d_out = [1,test_res];          % test_size is a random number

        case 'scalar'
            d_out = [1,1];

        otherwise
            error('Something went wrong in SOLUTION class method check_fcn_handle.');

    end

end



%% Do the check

% Check if function handle has the correct number of input arguments
if ~(nargin(fcn_handle)==1)     
    GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' has ',num2str(nargin(fcn_handle)),' arguments, but is only allowed to have 1.');


% function handle has the correct number of input arguments
else                            

    % Try to evaluate the function and check the sizes
    try                         

        % Evaluate the function with the test array ones(d_in)
        tmp = fcn_handle(ones(d_in));
        d_fcn = size(tmp);

        % Do the check for solget -> solget needs an extra check since the output size is variable and we only make sure that the size of every output dimension does not exceed the size of its corresponding input dimension
        if isfield(fcn_handle_options,'eval')
            if strcmpi(DYN.sol_type,'quasiperiodic') && (numel(d_fcn)==2) && (numel(d_in)==3)  % This is the case for quasi-periodic hypertime evaluations when output is not a 3D array
                d_fcn = [d_fcn, 1];                                                            % We need to add the 1 since Matlab automatically removes all array dimensions of size 1 "at the back"
            end
            if numel(d_fcn) > size(d_in,2)          % numel(d_fcn) <= size(d_in,2) is okay for solget
                GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns a ', num2str(numel(d_fcn)),'-dimensional array for a ', num2str(numel(d_in)),'-dimensional array input. However, the output should be a ', num2str(size(d_out,2)),'-dimensional array.');
            else
                for i = 1:numel(d_in)
                    if d_fcn(i) > d_in(i)
                        GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns for an exemplary input of size [', num2str(d_in),'] an output of size [', num2str(d_fcn),']. However, the maximum output size is [', num2str(d_in),'].');
                        break
                    end
                end
            end

        % Do the check for all other cases including contplot and solplot -> Output of function handle has a specific size, which we know. There is either exactly one possible size or there are dim possible sizes
        else    
            if strcmpi(DYN.sol_type,'quasiperiodic') && (numel(d_fcn)==2) && all(d_fcn == [test_res test_res])  % This is the case when 1 state variable is requested for hypertime plot of a quasi-periodic solution
                d_fcn = [d_fcn, 1];                                                                             % We need to add the 1 since Matlab automatically removes all array dimensions of size 1 "at the back"
            end
            if numel(d_fcn) ~= size(d_out,2)        % OLD: ~(numel(d_fcn)==numel(d_out)) -> can't be used anymore due to case 'matrix_and_vector'
                GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns a ', num2str(numel(d_fcn)),'-dimensional array for a ', num2str(numel(d_in)),'-dimensional array input. However, the output should be a ', num2str(size(d_out,2)),'-dimensional array.');
            elseif ~ismember(d_fcn,d_out,'rows')    % OLD: ~all(d_fcn==d_out) -> can't be used anymore due to case 'matrix_and_vector'
                if ~isvector(d_out)                 % This is for the case when multiple output sizes are possible (see e.g. 'matrix_or_vector')
                    GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns for an exemplary input of size [', num2str(d_in),'] an output of size [', num2str(d_fcn),']. However, the output size can only range from [', num2str(d_out(1,:)),'] to [',num2str(d_out(end,:)),'].');
                else
                    GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns for an exemplary input of size [', num2str(d_in),'] an output of size [', num2str(d_fcn),']. However, the output should be of size [', num2str(d_out),'].');
                end
            end

        end


    % Throw error message if function evaluation failed - this is (among other things) the case when the user requested array elements that do not exist (e.g. the user requested the third state variable via @(z) z(:,3), but the system has only two state variables)
    catch                       

        GC.error_msg{1,end+1} = append('Something with your provided function handle ',name_handle,' is wrong. It should have a structure similar to:  Input size: [', num2str(d_in),']  -->  Output size: [', num2str(d_out(1,:)),']');
        GC.error_msg{1,end+1} = 'Maybe you requested some array elements and/or dimensions that do not exist. Please be aware of the structure of your function handle, of the resolution (default value or your specified value) and of the system''s dimension.';
    
    end

end


end