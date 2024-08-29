% Function solget returns a solution at specific solution curve points 
%   [s_out,varargout] = solget(obj,DYN,options)
%
% Input Arguments:
% @obj:                 Object of Solution class
% @DYN:                 Object of DynamicalSystem class
% @options:             struct contain options for calculating solution
%
%                       Mandatory fields:
%                       'eval':     -->'function_handle:' User supplied
%                                       function handle. Example: options('eval',@(z)z(:,1))
%                                   --> 'euclidean': gives back the
%                                       euclidean (L2) along the states of an output vector
%                                       for the specified index.
%                                       Example: periodic solution
%                                       options('eval','euclidean') for an
%                                       200 x 3 x 10 (200 points per limit cycle, 3 states, 10 conti points)
%                                       a 200 x 1 x 10 vector
%                                   --> 'all': is not altering the output 
%
%                       'space':    --> 'solution': gives back the solution
%                                       vector(s). Example: s might be
%                                       shooting starting points or fourier
%                                       coefficients vectors for different
%                                       mu values
%                                   --> 'time': gives back the solution
%                                       over time. Example: Shooting:
%                                       Solution over time is (re)computed.
%                                       Fourier-Galerkin: Fourier
%                                       coefficients are evaluated for a
%                                       time vector
%                                   --> 'hypertime': gives back the
%                                       solution over hypertime. Periodic:
%                                       2D time trajectory over torus
%                                       coordinate \theta = [0,2\pi]
%                                       Quasi-periodic: 3D plane over
%                                       torus coordinate \theta_1,2 =
%                                       [0,2\pi]x [0, 2\pi]
%                                   --> 'frequency': gives back the FFT of
%                                       the hypertime solution (1D or nD FFT)
%
%                       Optional fields:
%                       'resolution'    1x1 or 1xnumb_base_frq double array
%                                       defining the number of points in
%                                       time or hypertime at which the
%                                       solution is evaluated. This is the
%                                       total number of points (see 'interval')
%                                       Default value: 200
%                                       Ignored for options.space = 'solution'  
%
%                       'interval'      1x2 interval: Defines the start and
%                                       end point of the evaluation time
%                                       interval 
%                                       Ignored for options.space =
%                                       'solution', 'hypertime'
%                                       Default value: [0,2\pi/frq]
%
%                       'index'         1xn array indiciating the indices
%                                       of the solution curve, where the
%                                       solution is evaluated.
%                                       Default: All solutions are given
%                                       back
%                       'mu'            1xn array indiciaing the mu values
%                                       of the solution curve, where the
%                                       solution is evaluated. If the exact
%                                       value is not present, the closest
%                                       value is used. If these leads to
%                                       doubled values, they are ignored.
%                                       Mu might lead to ambuguities in
%                                       case of overhanging curves
%
%
% Output Arguments: 
% @output:              Struct containing evaluated solution data, domain data, mu-values and options with the fields:
%                        - output.solution_eval: evaluated solution data
%                        - second field depends on solution space:
%                           * options.space = 'time'      --> output.time: time evaluation points
%                           * options.space = 'hypertime' --> output.hypertime: hypertime evaluation points
%                           * options.space = 'frequency' --> output.frequency: angular frequency values
%                        - output.mu: values of continuation parameter mu
%                        - output.options: options struct (can be equal to input options struct, but further fields might have been added in solget)
%
%
% Example:  options = costaropts('xaxis',@(z)z(:,1),'space','time','resolution',300,'interval',[0,2],'mu',[0.5,1,1.5]);
%           output = S.solget(DYN,options);


function output = solget(obj,DYN,options)
    
    %% Check the options structure: !!! BE CAREFUL - FUNCTION MAY CHANGE OPTIONS STRUCTURE AND SOLUTION OBJECT
    if ~isfield(options,'axes_values_old')              % When this field is present, contplot is called from the plot_contplot method during a continuation. The two following methods can be skipped in this case since the developers make sure that the options struct is fine
        options = obj.solget_gatekeeper(DYN,options);   % Check the input options structure
        options = obj.solget_up_index(DYN,options);     % Updates options.index: options.mu (if present) is replaced by options.index in this function. 
    end

    % The index is now unique, if mu leads to doubled values. options.index
    % can be used without caution for doubled indices or indices not
    % occurring in the SOLUTION subclass object

    % Set resolution default, if not defined
    if ~isfield(options,'resolution')
        if strcmpi(DYN.sol_type,'periodic')
            options.resolution = 200;
        elseif strcmpi(DYN.sol_type,'quasiperiodic')
            options.resolution = 50;
        end
    end


    %% Get the solutions and do the calculations (if necessary)

    switch options.space

        case 'time'
            [s,mu,time] = obj.evalsol_time(DYN,options);
            output.time = time;

        case 'hypertime'
            [s,mu,hypertime] = obj.evalsol_hypertime(DYN,options);
            output.hypertime = hypertime;

        case 'frequency'
            [s,alpha,mu,frequency] = obj.evalsol_frequency(DYN,options);
            output.frequency = frequency;
            output.angle = alpha;

    end


    % Get the correct array dimension of s (ATTENTION: ndims(s) returns 2 for matrices, vectors and scalars!)
    if iscolumn(s)                  % Column vector or scalar. It is important to check for column vector because EQ solutions can be a scalar and s therefore a [1x2] array (= vector) during continuation. In that case, s_array_dim has to be 2 because line 173 would not work otherwise
        s_array_dim = 1;            % Actually, 1 is not correct for scalar, but code below does not distinguish between vectors and scalars
    else
        s_array_dim = ndims(s);     % In all other cases, ndims(s) fits to the code below
    end


    % Is only one index / solution requested or is more than one index / solution requested ?
    idx_assist = double(~isscalar(options.index));      % If only one index / solution is requested, the second, third or fourth dimension of the output functions is only one / not present. For that case, idx_assist = 0.
                                                        % If more than one index / solution is requested, idx_assist = 1 and the correct case in the if-statement is chosen.


    % Apply the requested function to the data
    if isa(options.eval,'function_handle')

        if s_array_dim == (1+idx_assist)                % EQ 'hypertime'
            s_out = zeros([numel(options.eval(s(:,1))),numel(options.index)]);          % Initialise s_out
            for k = 1:numel(options.index)
                s_out(:,k) = options.eval(s(:,k));
            end
        elseif s_array_dim == (2+idx_assist)            % PS 'hypertime', 'time' and 'frequency'
            s_out = zeros([size(options.eval(s(:,:,1))),numel(options.index)]);         % Initialise s_out
            for k = 1:numel(options.index)
                s_out(:,:,k) = options.eval(s(:,:,k));
            end
        elseif s_array_dim == (3+idx_assist)            % QPS 'hypertime'
            s_out = zeros([size(options.eval(s(:,:,:,1)),[1 2]),size(options.eval(s(:,:,:,1)),3),numel(options.index)]); % Initialise s_out (third array dimension needs to be queried separately since it would be neglected by size() if it was equal to 1)
            for k = 1:numel(options.index)
               s_out(:,:,:,k) = options.eval(s(:,:,:,k));
            end
        end

    elseif ischar(options.eval)

        switch options.eval

            case 'euclidean'

               s_out = vecnorm(s,2,s_array_dim-idx_assist);  % s_array_dim-idx_assist is important since it assures that the norm is calculated within the correct dimension of s

            case 'all'    

               s_out = s;
               
        end

    end

    
    %% Output

    output.solution_eval = s_out;
    output.mu = mu;
    output.options = options;


end