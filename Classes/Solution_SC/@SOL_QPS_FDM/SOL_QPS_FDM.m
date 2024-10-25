% The subclass SOL_QPS_FDM inherits from Solution and adapts the class for periodic specific solutions approximated by finite-differences
% ATTENTION: This subclass must be adapted when error control for finite differences is implemented

classdef SOL_QPS_FDM < Solution

    properties

        freq            % Property to save the frequency
            
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    methods 
         
        %Abstract Superclass Method: Must be defined
        %
        % @y1:               solution curve point vector
        % @J1:               Jacobian Matrix w.r.t. solution space
        % @newton_flag:      newton_flag (duh) of Newton corrector
        % @DYN:              DynamicalSystem class object
        % @AM:               AM_QPS_FDM class object
        % @varargin:         Only filled if stability or/and error_control is activated   
        % @varargin{1,1}:    spectral error
        % @varargin{1,2}:    Floquet multipliers
        % @varargin{1,3}:    number of unstable Floquet multipliers
        % @varargin{1,4}:    largest absolute value of Floquet multiplier
        function IF_arch_init_data(obj,y1,J1,newton_flag,DYN,AM,varargin)       % Interface method for archiving the data of the initial solution

            obj.y0                  = y1;                                       
            obj.s(:,1)              = y1(1:(end-1-DYN.n_auto),1);               % Solution vector of approximation method without autonomous frequency
            obj.mu(1,1)             = y1(end,1);                                % Continuation parameter 
            obj.J{1,1}              = sparse(J1);                               % Jacobian matrix
            obj.dy(:,1)             = NaN(size(J1,1),1);                        % Initialised. Gets correctly filled by IF_arch_data
            obj.newton_flag(1,1)    = newton_flag;                              % Exit-flag of corrector (fsolve)
            obj.arclength(1,1)      = 0;                                        % Set arclength of first curve point to zero
            
            if DYN.n_auto == 0
                obj.freq(:,1) = reshape(DYN.non_auto_freq(y1(end,1)),2,1);      % Frequencies if system is non-autonomous
            elseif DYN.n_auto == 1
                obj.freq(:,1) = [DYN.non_auto_freq(y1(end,1)); y1(end-1,1)];    % Frequencies if system is partly autonomous
            elseif DYN.n_auto == 2
                obj.freq(:,1) = y1(end-2:end-1,1);                              % Frequencies if system is fully autonomous
            end

            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,1)    = varargin{1,1}{1,2};
                obj.vectors(:,:,1)      = varargin{1,1}{1,5};
                obj.n_unstable(1,1)     = varargin{1,1}{1,3};
            end

        end
        
        %Abstract Superclass Method: Must be defined
        %
        %@CON:  Continuation class object
        %@DYN:  DynamicalSystem class object
        %@AM:   AM_QPS_FDM class object
        function IF_arch_data(obj,CON,DYN,AM)                                       % Interface method for archiving the data of the continuation

            obj.s(:,end+1)              = CON.p_y1(1:(end-1-DYN.n_auto),1);         % Solution vector of approximation method without autonomous frequency
            obj.mu(1,end+1)             = CON.p_y1(end,1);                          % Continuation parameter
            obj.J{1,end+1}              = sparse(CON.p_J1);                         % Jacobian matrix
            obj.dy(:,end:end+1)         = [CON.dy0, NaN(size(CON.dy0))];            % Direction vector of the predictor
            obj.newton_flag(1,end+1)    = CON.p_newton_flag;                        % Exit-flag of corrector (fsolve)
            obj.step_width(1,end+1)     = CON.step_width;                           % Step width
            obj.arclength(1,end+1)      = CON.p_arcl_1;                             % Arclength

            if DYN.n_auto == 0
                obj.freq(:,end+1) = reshape(DYN.non_auto_freq(CON.p_y1(end,1)),2,1);            % Frequencies if system is non-autonomous
            elseif DYN.n_auto == 1
                obj.freq(:,end+1) = [DYN.non_auto_freq(CON.p_y1(end,1)); CON.p_y1(end-1,1)];    % Frequencies if system is partly autonomous
            elseif DYN.n_auto == 2
                obj.freq(:,end+1) = CON.p_y1(end-2:end-1,1);                                    % Frequencies if system is fully autonomous
            end

            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,end+1)    = CON.p_multipliers;                    % 
                obj.vectors(:,:,end+1)      = CON.p_vectors;                        % 
                obj.n_unstable(1,end+1)     = CON.p_n_unstable_1;                   % Number of unstable multipliers
            end

        end

        %Abstract Superclass Method: Must be defined
        %
        %@CON:  Continuation class object
        %@DYN:  DynamicalSystem class object
        %@AM:   AM_QPS_FDM class object
        function IF_arch_bfp_data(obj,CON,DYN,AM,ST)    % NOT WORKING CURRENTLY: Interface method for archiving the data of an iterated bifurcation point

            % obj.s(:,end+1)          = CON.p_y_bfp(1:(end-1-DYN.n_auto),1);                      % Solution vector of approximation method
            % obj.mu(1,end+1)         = CON.p_y_bfp(end,1);                                       % Continuation parameter
            % obj.J(:,:,end+1)        = CON.p_J_bfp;                                              % Jacobian matrix
            % obj.dy(:,end:end+1)     = [NaN(size(CON.dy0)), NaN(size(CON.dy0))];                 % Direction vector of the predictor
            % obj.newton_flag(1,end+1)= NaN;                                                      % Exit-flag is unknown as it is not saved as a Stability class property 
            % obj.step_width(1,end+1) = CON.step_width;                                           % Step width
            % obj.arclength(1,end+1)  = CON.p_arclength_bfp;                                      % Arc length
            % 
            % if DYN.n_auto == 0
            %     obj.freq(:,end+1) = reshape(DYN.non_auto_freq(CON.p_y_bfp(end,1)),2,1);             % Frequencies if system is non-autonomous
            % elseif DYN.n_auto == 1
            %     obj.freq(:,end+1) = [DYN.non_auto_freq(CON.p_y_bfp(end,1)); CON.p_y_bfp(end-1,1)];  % Frequencies if system is partly autonomous
            % elseif DYN.n_auto == 2
            %     obj.freq(:,end+1) = CON.p_y_bfp(end-2:end-1,1);                                     % Frequencies if system is fully autonomous
            % end
            % 
            % obj.multipliers(:,end+1)    = CON.p_multipliers;                % 
            % obj.vectors(:,:,end+1)      = CON.p_vectors_bfp;                % 
            % obj.n_unstable(1,end+1)     = obj.n_unstable(1,end);            % Indiacting number of unstable multipliers. Definition: The number in the point is equal to the number before the bfp
            % 
            % % Fill the table for the bifurcations 
            % [label,msg] = ST.identify_bifurcation();
            % obj.bifurcation = [obj.bifurcation;{label,numel(obj.mu),msg}];
        
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
        
    methods(Access = protected)         % Can only be called by within class and superclass objects

        %Abstract Superclass Methods: Must be defined with the exact same passed and returned arguments
        
        %Postprocessing
        [s_time,mu,time]                        = evalsol_time(obj,DYN,options);                %Function gives back the solution in time domain
        [s_hypertime,mu,hypertimes]             = evalsol_hypertime(obj,DYN,options);           %Function gives back the solution in hypertime domain
        [s_amplitude,s_angle,mu,frequency]      = evalsol_frequency(obj,DYN,options);           %Function gives back the solution in frequency domain

    end

end