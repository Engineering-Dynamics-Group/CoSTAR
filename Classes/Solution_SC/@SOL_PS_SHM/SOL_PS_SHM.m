%The subclass SOL_PS_SHM inherits from Solution and adapts the class for periodic solutions approximated by using the shooting method

classdef SOL_PS_SHM < Solution

    properties

        freq                                    % storing frequency values (they are either non-autonomous or autonomous for periodic solutions)
        solver_function function_handle
        odeOpts

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        % Constructor
        function obj = SOL_PS_SHM(AM)
           obj.solver_function = AM.solver_function;
           obj.odeOpts = AM.odeOpts;
        end
        
        % Interface method for archiving the data of the initial solution
        function IF_arch_init_data(obj,y1,J1,newton_flag,DYN,AM,varargin)                   
            
            obj.y0                     = y1;
            obj.s(:,1)                 = y1(1:(end-1-DYN.n_auto),1);                    % Approximation method vector 
            obj.mu(1,1)                = y1(end,1);                                     % Continuation parameter 
            obj.J{1,1}                 = sparse(J1);                                    % Jacobian matrix
            obj.dy(:,1)                = NaN(size(J1,1),1);                             % Initialised. Gets correctly filled by IF_arch_data
            obj.newton_flag(1,1)       = newton_flag;
            obj.arclength(1,1)         = 0;                                             % Set arc length of first curve point to zero
            
            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,1)    = varargin{1,1}{1,2};
                obj.vectors(:,:,1)      = varargin{1,1}{1,5};
                obj.n_unstable(1,1)     = varargin{1,1}{1,3}; 
                obj.stability_flag(1,1) = varargin{1,1}{1,4};
            end

            if DYN.n_auto == 0                                                          % Solution is periodic - if this is true: non-autonomous
                obj.freq(1,1) = DYN.non_auto_freq(y1(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,1) = y1(end-1,1);
            end

        end
        
        % Interface method for archiving the data of the continuation
        function IF_arch_data(obj,CON,DYN,AM)                                                  
            
            obj.s(:,end+1)          = CON.p_y1(1:(end-1-DYN.n_auto),1);                 % Approximation method vector 
            obj.mu(1,end+1)         = CON.p_y1(end,1);                                  % Continuation parameter 
            obj.J{1,end+1}          = sparse(CON.p_J1);                                 % Jacobian matrix
            obj.dy(:,end:end+1)     = [CON.dy0, NaN(size(CON.dy0))];                    % Direction vector of the predictor
            obj.newton_flag(1,end+1)= CON.p_newton_flag;                                % Exit-flag of corrector (fsolve)
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arcl_1;
            
            if DYN.n_auto == 0                                                          % Solution is periodic - if this is true: non-autonomous
                obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_y1(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,end+1) = CON.p_y1(end-1,1);
            end

            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,end+1)    = CON.p_multipliers;                        % Floquet Multipliers
                obj.vectors(:,:,end+1)      = CON.p_vectors;                            % Eigenvectors corresponding to Floquet Multipliers
                obj.n_unstable(1,end+1)     = CON.p_n_unstable_1;                       % Indiacting number of unstable multipliers
                obj.stability_flag(1,end+1) = CON.p_stability_flag;                     %Exitflag of stability computation
            end

        end

        % Interface method for archiving the data of a bifurcation point
        function IF_arch_bfp_data(obj,CON,DYN,AM,ST)

            obj.s(:,end+1)          = CON.p_y_bfp(1:(end-1-DYN.n_auto),1);              % Approximation method vector
            obj.mu(1,end+1)         = CON.p_y_bfp(end,1);                               % Continuation parameter
            obj.J{1,end+1}          = sparse(CON.p_J_bfp);                              % Jacobian matrix
            obj.dy(:,end:end+1)     = [NaN(size(CON.dy0)), NaN(size(CON.dy0))];         % Direction vector of the predictor
            obj.newton_flag(1,end+1)= CON.p_newton_flag_bfp;                            % Exit-flag of corrector (fsolve)
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arclength_bfp;

            if DYN.n_auto == 0                                                          % Solution is periodic - if this is true: non-autonomous
                obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_y_bfp(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,end+1) = CON.p_y_bfp(end-1,1);
            end

            obj.multipliers(:,end+1)    = CON.p_multipliers_bfp;                        % Floquet Multipliers
            obj.vectors(:,:,end+1)      = CON.p_vectors_bfp;                            % Eigenvectors corresponding to Floquet Multipliers
            obj.n_unstable(1,end+1)     = obj.n_unstable(1,end);                        % Indiacting number of unstable multipliers.
                                                                                        % Definition: The number in the point is equal to the number before the bfp
            obj.stability_flag(1,end+1) = CON.p_stability_flag;                         % Exitflag of stability computation
                                                                                        % Fill the table for the bifurcations 
            [label,msg] = ST.identify_bifurcation();
            obj.bifurcation = [obj.bifurcation;{label,numel(obj.mu),msg}];

        end

    end

    %%%%%%%%%%%%%%%
        
    methods(Access = protected)     % Can only be called from within class and superclass objects

        % Postprocessing
        [s_time,mu,time]                        = evalsol_time(obj,DYN,options);            % Function returning the solution in time domain
        [s_hypertime,mu,hypertimes]             = evalsol_hypertime(obj,DYN,options);       % Function returning the solution in hypertime domain
        [s_amplitude,s_angle,mu,frequency]      = evalsol_frequency(obj,DYN,options);       % Function returning the solution in frequency domain

    end

end