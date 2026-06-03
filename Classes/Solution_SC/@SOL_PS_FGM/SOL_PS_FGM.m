%The subclass SOL_PS_FGM inherits from Solution and adapts the class for
%periodic  specific solutions computed with the Fourier-Galerkin method

classdef SOL_PS_FGM < Solution

    properties

        freq            %Property to save the frequency
        hmatrix         %Property to save the higher-harmonics matrix for every continuation step
        n_hh            %Number of higher harmonics
        n_fft           %Number of FFT evaluation points
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    methods 

        %Abstract Superclass Method: Must be defined
        %Method for storing the initial solution into the solution class object
        %
        %@y1:               solution curve point vector
        %@J1:               Jacobian Matrix w.r.t. solution space
        %@newton_flag:      newton_flag (duh) of Newton corrector
        %@DYN:              DynamicalSystem class object
        %@AM:               AM_PS_FGM class object
        %@varargin:         Only filled if stability or/and error_control is activated   
        %@varargin{1,1}:    spectral error
        %@varargin{1,2}:    Floquet multipliers
        %@varargin{1,3}:    number of unstable Floquet multipliers
        %@varargin{1,4}:    largest absolute value of Floquet multiplier
        function IF_arch_init_data(obj,y1,J1,newton_flag,DYN,AM,varargin)                               %Interface method for archiving the data of the initial solution
                
               %This method gets called by the intial solution function

               %Store here the intial solution properties of the superclass and
               %your newly defined subclass properties
  
            obj.y0                     = y1;
            obj.s{1,1}                 = y1(1:(end-1-DYN.n_auto),1);                        %For equilibria, there is no autonomous frequency in the solution vector
            obj.mu(1,1)                = y1(end,1);                             
            obj.J{1,1}                 = J1;
            obj.dy{1,1}                = NaN(size(y1));                                     %Initialised. Gets correctly filled by IF_arch_data
            obj.newton_flag(1,1)       = newton_flag;
            obj.arclength(1,1)         = 0;                                                 %Set arclength of first curve point to zero
            obj.hmatrix{1,1}           = AM.hmatrix;
            obj.n_hh(1,1)              = numel(AM.hmatrix);
            obj.n_fft(1,1)             = AM.n_fft;

            if strcmpi(AM.error_control,'on'); obj.error(1,1)      = varargin{1,1}{1,1}; end
            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,1)    = varargin{1,1}{1,2};
                obj.vectors(:,:,1)      = varargin{1,1}{1,5};
                obj.n_unstable(1,1)     = varargin{1,1}{1,3}; 
                obj.stability_flag(1,1) = varargin{1,1}{1,4};
            end
    


            if DYN.n_auto == 0                                                              %Solution is periodic - if this is true: non-autonomous
                obj.freq(1,1) = DYN.non_auto_freq(y1(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,1) = y1(end-1,1);
            end


        end
        
        %Abstract Superclass Method: Must be defined
        %Method for storing the current solution point into the solution class object
        %
        %@CON:  Continuation class object
        %@DYN:  DynamicalSystem class object
        %@AM:   AM_PS_FGM class object

        function IF_arch_data(obj,CON,DYN,AM)                                                  %Interface method for archiving the data of the continuation
     
               %This method gets called by the m_continuation method

               %Store here the solution properties of the superclass and
               %your newly defined subclass properties
  
               %Example:

            obj.s{1,end+1}          = CON.p_y1(1:(end-1-DYN.n_auto),1);                         %approximation method vector 
            obj.mu(1,end+1)         = CON.p_y1(end,1);                                          %continuation parameter 
            obj.J{1,end+1}          = CON.p_J1;                                                 %Jacobian matrix
            obj.dy{1,end}           = CON.dy0;                                                  %direction vector of the predictor (this one belongs to the previous point)
            obj.dy{1,end+1}         = NaN(size(CON.dy0));                                       %expand dy by a NaN vector which can be filled with the direction vector in the next loop
            obj.newton_flag(1,end+1)= CON.p_newton_flag;
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arcl_1;
            obj.hmatrix{1,end+1}    = AM.hmatrix;
            if strcmpi(AM.error_control,'on')
                obj.error(1,end+1)  = CON.p_error;
            end
            obj.n_hh(1,end+1)       = numel(AM.hmatrix);
            obj.n_fft(1,end+1)      = AM.n_fft;

            if DYN.n_auto == 0                                                              %Solution is periodic - if this is true: non-autonomous
                obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_y1(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,end+1) = CON.p_y1(end-1,1);
            end

            if strcmpi(DYN.stability,'on')
                obj.multipliers(:,end+1)    = CON.p_multipliers;                        %Floquet Multipliers
                obj.vectors(:,:,end+1)      = CON.p_vectors;                            %Eigenvectors corresponding to Floquet Multipliers
                obj.n_unstable(1,end+1)     = CON.p_n_unstable_1;                       %Indiacting number of unstable multipliers
                obj.stability_flag(1,end+1) = CON.p_stability_flag;                     %Exitflag of stability computation
            end

     

        end

        %Abstract Superclass Method: Must be defined
        %Method for storing an iterated bifurcation point in the solution class object
        %
        %@CON:  Continuation class object
        %@DYN:  DynamicalSystem class object
        %@AM:   AM_PS_FGM class object

        function IF_arch_bfp_data(obj,CON,DYN,AM,ST)

            obj.s{1,end+1}          = CON.p_y_bfp(1:(end-1-DYN.n_auto),1);                      %approximation method vector
            obj.mu(1,end+1)         = CON.p_y_bfp(end,1);                                       %continuation parameter
            obj.J{1,end+1}          = CON.p_J_bfp;                                              %Jacobian matrix
            obj.dy{1,end}           = NaN(size(CON.dy0));                                       %direction vector of the predictor (this one belongs to the previous point)
            obj.dy{1,end+1}         = NaN(size(CON.dy0));                                       %expand dy by a NaN vector which can be filled with the direction vector in the next loop
            obj.newton_flag(1,end+1)= CON.p_newton_flag_bfp;                                    %Exit-flag of corrector (fsolve)
            obj.step_width(1,end+1) = CON.step_width;
            obj.arclength(1,end+1)  = CON.p_arclength_bfp;
            obj.hmatrix{1,end+1}    = AM.hmatrix;                                               %This is the same as the one of the newest point (due to a update_sol_dim statement in the ST.update_curve_container method)                                       
            obj.n_hh(1,end+1)       = numel(AM.hmatrix);
            obj.n_fft(1,end+1)      = AM.n_fft;
            if strcmpi(AM.error_control,'on')
                obj.error(1,end+1)  = CON.p_error_bfp;                                          %error at the bifurcation point
            end

            if DYN.n_auto == 0                                                                   %Solution is periodic - if this is true: non-autonomous
                obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_y_bfp(end,1));
            elseif DYN.n_auto == 1
                obj.freq(1,end+1) = CON.p_y_bfp(end-1,1);
            end

            obj.multipliers(:,end+1)    = CON.p_multipliers_bfp;                    %Floquet Multipliers
            obj.vectors(:,:,end+1)      = CON.p_vectors_bfp;                        %Eigenvectors corresponding to Floquet Multipliers
            obj.n_unstable(1,end+1)     = obj.n_unstable(1,end);                    %Indiacting number of unstable multipliers. Definition: The number in the point is equal to the number before the bfp 
            obj.stability_flag(1,end+1) = CON.p_stability_flag;                     %Exitflag of stability computation

            %Fill the table for the bifurcations 
            [label,msg] = ST.identify_bifurcation();
            obj.bifurcation = [obj.bifurcation;{label,numel(obj.mu),msg}];

        end

    end

        
    methods(Access = protected) %Can only be called by within class and superclass objects

        %Abstract Superclass Methods: Must be defined with the exact same
        %passed and returned arguments
        
        %Postprocessing
        [s_time,mu,time]                        = evalsol_time(obj,DYN,options);                %Function gives back the solution in time domain
        [s_hypertime,mu,hypertimes]             = evalsol_hypertime(obj,DYN,options);           %Function gives back the solution in hypertime domain
        [s_amplitude,s_angle,mu,frequency]      = evalsol_frequency(obj,DYN,options);           %Function gives back the solution in frequency domain

    end

end