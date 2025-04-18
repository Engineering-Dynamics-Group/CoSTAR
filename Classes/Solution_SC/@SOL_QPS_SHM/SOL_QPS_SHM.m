%The subclass SOL_QPS_SHM inherits from Solution and adapts the class for
%quasi-periodic specific solutions calculated by a shooting algorithm

classdef SOL_QPS_SHM < Solution
    
    properties
        
        freq                                                                                % storing frequency values (they are either non-autonomous or autonomous for periodic solutions)
        solver_function function_handle                                                     % Ode solver function (i.e. ode45, ode15s)
        odeOpts                                                                             % Options for ode solver
        Ik                                                                                  % Integration invterval
        phi                                                                                 % Phase values for each characteristic
        n_char                                                                              % Number of characteristics
        n                                                                                   % Number of state-space variables
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function obj = SOL_QPS_SHM(AM)
            obj.solver_function = AM.solver_function;                                       % Fetch solver-function (i.e. ode45)
            obj.odeOpts = AM.odeOpts;                                                       % Fetch ode-solver options
            obj.n_char = AM.n_char;                                                         % Fetch number of characteristics
            obj.n = AM.n;                                                                   % Fetch number of state-space variables
        end
        
        function IF_arch_init_data(obj,y1,J1,Ik,newton_flag,phi,DYN,varargin)               % Interface method for archiving the data of the initial solution
            
            obj.y0                     = y1;
            obj.s(:,1)                 = y1(1:(end-1-DYN.n_auto),1);                        % Set only the solutions points not the autonous frequencies and bifurcation parameter
            obj.mu(1,1)                = y1(end,1);                                         % Bifurcation parameter
            obj.J{1,1}                 = sparse(J1);                                        % Jacobian
            obj.dy(:,1)                = NaN(size(J1,1),1);                                 % Initialised. Gets correctly filled by IF_arch_data
            obj.newton_flag(1,1)       = newton_flag;                                       % Exitflag of fsolve
            obj.arclength(1,1)         = 0;                                                 % Set arclength of first curve point to zero
            obj.Ik(1,:)                = Ik;                                                % Set integration time
            obj.phi                    = phi;                                               
            
            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,1)    = varargin{1,1}{1,2};
                obj.vectors(:,:,1)      = varargin{1,1}{1,5};
                obj.n_unstable(1,1)     = varargin{1,1}{1,3}; 
                obj.stability_flag(1,1) = varargin{1,1}{1,4};
            end

            if DYN.n_auto == 0                                      
                obj.freq(:,end+1) = DYN.non_auto_freq(y1(end,1));                           % Non-autonomous case both frequencies by function DYN.non_auto_freq
            elseif DYN.n_auto == 1
                obj.freq(:,end+1) = [DYN.non_auto_freq(y1(end,1));y1(end-1,1)];             % Mixed case non-autonomous frequency by function DYN.non_auto_freq, autonomous frequency by solution vector
            elseif DYN.n_auto == 2
                obj.freq(:,end+1) = [y1(end-2,1);y1(end-1,1)];                              % Full-autonomous case both autonomous frequency by solution vector
            end
            
        end
        
        function IF_arch_data(obj,CON,DYN,AM)                                               %Interface method for archiving the data of the continuation
            
            obj.s(:,end+1)          = CON.p_y1(1:(end-1-DYN.n_auto),1);                     % solution method vector
            obj.mu(1,end+1)         = CON.p_y1(end,1);                                      % continuation parameter
            obj.J{1,end+1}          = sparse(CON.p_J1);                                     % Jacobian matrix
            obj.dy(:,end:end+1)     = [CON.dy0, NaN(size(CON.dy0))];                        % Direction vector of the predictor
            obj.newton_flag(1,end+1)= CON.p_newton_flag;                                    % Exitflag of fsolve
            obj.step_width(1,end+1) = CON.step_width;                                       % Current step width
            obj.arclength(1,end+1)  = CON.p_arcl_1;                                         % Current value of arclength
            
            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,end+1)    = CON.p_multipliers;                            %
                obj.vectors(:,:,end+1)      = CON.p_vectors;                                %
                obj.n_unstable(1,end+1)     = CON.p_n_unstable_1;                           % Indiacting number of unstable multipliers
                obj.stability_flag(1,end+1) = CON.p_stability_flag;                         % Exitflag of stability computation
            end

            if DYN.n_auto == 0                                         
                obj.freq(:,end+1) = DYN.non_auto_freq(CON.p_y1(end,1));                     % Non-autonomous case both frequencies by function DYN.non_auto_freq
            elseif DYN.n_auto == 1
                obj.freq(:,end+1) = [DYN.non_auto_freq(CON.p_y1(end,1));CON.p_y1(end-1,1)]; % Mixed case non-autonomous frequency by function DYN.non_auto_freq, autonomous frequency by solution vector
            elseif DYN.n_auto == 2
                obj.freq(:,end+1) = [CON.p_y1(end-2,1);CON.p_y1(end-1,1)];                  % Full-autonomous case both autonomous frequency by solution vector
            end 
        end

        function IF_arch_bfp_data(obj,CON,DYN,AM,ST)

            obj.s(:,end+1)          = CON.p_y_bfp(1:(end-1-DYN.n_auto),1);                      %approximation method vector
            obj.mu(1,end+1)         = CON.p_y_bfp(end,1);                                       %continuation parameter
            obj.J{1,end+1}          = sparse(CON.p_J_bfp);                                      %Jacobian matrix
            obj.dy(:,end:end+1)     = [NaN(size(CON.dy0)), NaN(size(CON.dy0))];                 %Direction vector of the predictor
            obj.newton_flag(1,end+1)= CON.p_newton_flag_bfp;                                    %Exit-flag of corrector (fsolve);                                                      %Exit-flag is unknown as it is not saved as a Stability class property
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arclength_bfp;

            if DYN.n_auto == 0                                         
                obj.freq(:,end+1) = DYN.non_auto_freq(CON.p_y1(end,1));                     % Non-autonomous case both frequencies by function DYN.non_auto_freq
            elseif DYN.n_auto == 1
                obj.freq(:,end+1) = [DYN.non_auto_freq(CON.p_y1(end,1));CON.p_y1(end-1,1)]; % Mixed case non-autonomous frequency by function DYN.non_auto_freq, autonomous frequency by solution vector
            elseif DYN.n_auto == 2
                obj.freq(:,end+1) = [CON.p_y1(end-2,1);CON.p_y1(end-1,1)];                  % Full-autonomous case both autonomous frequency by solution vector
            end 

            obj.multipliers(:,end+1)    = CON.p_multipliers_bfp;                    % Ljapunov exponents
            obj.vectors(:,:,end+1)      = CON.p_vectors_bfp;                        % There are no vectors related to Ljapunov exponents, so this is empty
            obj.n_unstable(1,end+1)     = obj.n_unstable(1,end);                    % Indiacting number of unstable multipliers. Definition: The number in the point is equal to the number before the bfp 
            obj.stability_flag(1,end+1) = CON.p_stability_flag;                     % Exitflag of stability computation

            %Fill the table for the bifurcations 
            [label,msg] = ST.identify_bifurcation();
            obj.bifurcation = [obj.bifurcation;{label,numel(obj.mu),msg}];

        end
        
    end
    
    
    methods(Access = protected) %Can only be called by within class and superclass objects
        
        %Postprocessing
        [s_time,mu,time]                        = evalsol_time(obj,DYN,options);            % Function gives back the solution in time domain
        [s_hypertime,mu,hypertimes]             = evalsol_hypertime(obj,DYN,options);       % Function gives back the solution in hypertime domain
        [s_amplitude,s_angle,mu,frequency]      = evalsol_frequency(obj,DYN,options);       % Function gives back the solution in frequency domain
        
        f                                       = FcnWrapper_SOL_ODE2(obj,t,z,Fcn,PHI)      % Function wrapper for time integration of all characteristics
    end
    
end