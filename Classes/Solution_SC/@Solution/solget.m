% Function solget returns a solution at specific solution curve points 
%   [s_out,varargout] = solget(obj,DYN,options)
%
% Input Arguments
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
% Output Arguments: 
% @s_out:               Solution vector in the dimension [options.resolution,
%                       statespace dimension, number of curve points].
%                       statespace dimension depends on options.eval    
% @varargout:           varargout{1,1}: mu - array of continuation parameters, 
%                       where solution was truely evaluated
%                       varargout{1,2}: x - time/hypertime/frequency array. 
%                       for option.space = 'solution'. Nothing is given
%                       back.
%                       varargout{1,3}: options structure where options.mu
%                       might have been replaced by index.
%
%Example:   options = struct('xaxis',@(z)z(:,1),'space','time','resolution',300,'interval',[0,2],'mu',[0.5,1,1.5]);
%           [s,mu,time] = S.solget(DYN,options);


function [s_out,mu,x,options] = solget(obj,DYN,options)            % [s_out,varargout]  
    
    %% Check the options structure: !!! BE CAREFUL - FUNCTION MAY CHANGE OPTIONS STRUCTURE AND SOLUTION OBJECT
    if ~isfield(options,'axes_values_old')              %When this field is present, contplot is called from the plot_contplot method during a continuation. The two following methods can be skipped in this case since the developers make sure that the options struct is fine
        options = obj.solget_gatekeeper(DYN,options);   %Check the input options structure
        options = obj.solget_up_index(DYN,options);     %Updates options.index: options.mu (if present) is replaced by options.index in this function. 
    end

    % The index is now unique, if mu leads to doubled values. options.index
    % can be used without caution for doubled indices or indices not
    % occurring in the SOLUTION subclass object

    %Set resolution default, if not defined
    if ~isfield(options,'resolution')
        options.resolution = 200;
    end

    %% Get the solutions and do the calculations (if necessary)
    % varargout = cell(1,3);

    switch options.space

        case 'time'
            %This method must return array structures
            [s,mu,x]  = obj.evalsol_time(DYN,options);
            % varargout{1,2} = x;

        case 'hypertime' %time and trajectory store to the same thing... mainly relevant for 
            %This method must return array structures
            [s,mu,x]   = obj.evalsol_hypertime(DYN,options);
            % varargout{1,2} = x;

        case 'frequency'        
            %This method must return array structures: It operates on
            %hyper-time series.
            [s,mu,x]   = obj.evalsol_frequency(DYN,options);
            % varargout{1,2} = x;

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
            for k = 1:numel(options.index)
                s_out(:,k) = options.eval(s(:,k));
            end
        elseif s_array_dim == (2+idx_assist)            % PS 'hypertime', 'time' and 'frequency'
            for k = 1:numel(options.index)
                s_out(:,:,k) = options.eval(s(:,:,k));
            end
        elseif s_array_dim == (3+idx_assist)            % QPS 'hypertime'
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

    % varargout{1,1} = mu;
    % varargout{1,3} = options;


end