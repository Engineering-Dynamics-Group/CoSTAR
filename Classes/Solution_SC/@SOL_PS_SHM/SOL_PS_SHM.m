%The subclass SOL_PS_SHM inherits from Solution and adapts the class for
%periodic  specific solutions

classdef SOL_PS_SHM < Solution

    properties

            freq                                                            %storing frequency values (they are either non-autonomous or autonomous for periodic solutions)
           solver_function function_handle
           odeOpts
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    methods 
        function obj = SOL_PS_SHM(AM)
           obj.solver_function = AM.solver_function;
           obj.odeOpts = AM.odeOpts;
        end
        
        function IF_arch_init_data(obj,y1,J1,newton_flag,DYN,AM,varargin)                   %Interface method for archiving the data of the initial solution
            
            obj.y0                     = y1;
            obj.s(:,1)                 = y1(1:(end-1-DYN.n_auto),1);                        %For equilibria, there is no autonomous frequency in the solution vector
            obj.mu(1,1)                = y1(end,1);                             
            obj.J(:,:,1)               = J1;
            obj.dy(:,1)                = NaN(size(J1,1),1);                                 %Initialised. Gets correctly filled by IF_arch_data
            obj.newton_flag(1,1)       = newton_flag;
            obj.arclength(1,1)         = 0;                                                 %Set arclength of first curve point to zero
            
            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,1)    = varargin{1,1}{1,2};
                obj.vectors(:,:,1)      = varargin{1,1}{1,5};
                obj.n_unstable(1,1)     = varargin{1,1}{1,3}; 
            end


            if DYN.n_auto == 0                                                              %Solution is periodic - if this is true: non-autonomous
                obj.freq(1,1) = DYN.non_auto_freq(y1(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,1) = y1(end-1,1);
            end

        end
        
        function IF_arch_data(obj,CON,DYN,AM)                                                  %Interface method for archiving the data of the continuation
            
            obj.s(:,end+1)          = CON.p_y1(1:(end-1-DYN.n_auto),1);                       %approximation method vector 
            obj.mu(1,end+1)         = CON.p_y1(end,1);                                        %continuation parameter 
            obj.J(:,:,end+1)        = CON.p_J1;                                               %Jacobian matrix
            obj.dy(:,end:end+1)     = [CON.dy0, NaN(size(CON.dy0))];                          %Direction vector of the predictor
            obj.newton_flag(1,end+1)= CON.p_newton_flag;                                      %Exit-flag of corrector (fsolve)
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arcl_1;
            
            if DYN.n_auto == 0                                                              %Solution is periodic - if this is true: non-autonomous
                obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_y1(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,end+1) = CON.p_y1(end-1,1);
            end

            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,end+1)    = CON.p_multipliers;                        %Floquet Multipliers
                obj.vectors(:,:,end+1)      = CON.p_vectors;                            %Eigenvectors corresponding to Floquet Multipliers
                obj.n_unstable(1,end+1)     = CON.p_n_unstable_1;                       %Indiacting number of unstable multipliers
            end

        end



        function IF_arch_bfp_data(obj,CON,DYN,AM,ST)

            obj.s(:,end+1)          = CON.p_y_bfp(1:(end-1-DYN.n_auto),1);                      %approximation method vector
            obj.mu(1,end+1)         = CON.p_y_bfp(end,1);                                       %continuation parameter
            obj.J(:,:,end+1)        = CON.p_J_bfp;                                              %Jacobian matrix
            obj.dy(:,end:end+1)     = [NaN(size(CON.dy0)), NaN(size(CON.dy0))];                 %Direction vector of the predictor
            obj.newton_flag(1,end+1)= NaN;                                                      %Exit-flag is unknown as it is not saved as a Stability class property
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arclength_bfp;

            if DYN.n_auto == 0                                                              %Solution is periodic - if this is true: non-autonomous
                obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_y_bfp(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,end+1) = CON.p_y_bfp(end-1,1);
            end

            obj.multipliers(:,end+1)    = CON.p_multipliers_bfp;                    %Floquet Multipliers
            obj.vectors(:,:,end+1)      = CON.p_vectors_bfp;                        %Eigenvectors corresponding to Floquet Multipliers
            obj.n_unstable(1,end+1)     = obj.n_unstable(1,end);                    %Indiacting number of unstable multipliers. Definition: The number in the point is equal to the number before the bfp 
       
            %Fill the table for the bifurcations 
            [label,msg] = ST.identify_bifurcation();
            obj.bifurcation = [obj.bifurcation;{label,numel(obj.mu),msg}];

        end

    end



        
    methods(Access = protected) %Can only be called by within class and superclass objects

        %Postprocessing
        [s_time,mu,time]                        = evalsol_time(obj,DYN,options);                %Function gives back the solution in time domain
        [s_hypertime,mu,hypertimes]             = evalsol_hypertime(obj,DYN,options);           %Function gives back the solution in hypertime domain
        [s_amplitude,s_angle,mu,frequency]      = evalsol_frequency(obj,DYN,options);           %Function gives back the solution in frequency domain

    end

end