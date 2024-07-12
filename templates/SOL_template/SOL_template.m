%The subclass PS_Shoot_Sol inherits from Solution and adapts the class for
%periodic  specific solutions

classdef SOL_template < Solution

    properties

        %Most of the relevant solution properties are defined in the
        %Solution superclass. Define here only the needed extensions
        %
        %e.g. a property to save the frequency
            
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    methods 

        %Abstract Superclass Method: Must be defined
        function IF_arch_init_data(obj,y1,J1,newton_flag,DYN,AM,varargin)                               %Interface method for archiving the data of the initial solution
                
               %This method gets called by the intial solution function

               %Store here the intial solution properties of the superclass and
               %your newly defined subclass properties
  
               %Example:

%             obj.y0                     = y1;
%             obj.s(:,1)                 = y1(1:(end-1-DYN.n_auto),1);                        %For equilibria, there is no autonomous frequency in the solution vector
%             obj.mu(1,1)                = y1(end,1);                             
%             obj.J(:,:,1)               = J1;
%             obj.dy(:,1)                = zeros(size(J1,1),1);                               %Initialised. Gets correctly filled by IF_arch_data
%             obj.newton_flag(1,1)       = newton_flag;
%             obj.arclength(1,1)         = 0;                                                 %Set arclength of first curve point to zero
%            if strcmpi(AM.error_control,'on'); obj.error(1,1)      = varargin{1,1}{1,1}; end
%             if strcmpi(DYN.stability,'on')
%                 obj.multipliers(:,1)    = varargin{1,1}{1,2};
%                 obj.n_unstable(1,1)     = varargin{1,1}{1,3}; 
%             end

        end
        
        %Abstract Superclass Method: Must be defined
        function IF_arch_data(obj,CON,DYN,AM)                                                  %Interface method for archiving the data of the continuation
     
               %This method gets called by the m_continuation method

               %Store here the solution properties of the superclass and
               %your newly defined subclass properties
  
               %Example:

%             obj.s(:,end+1)          = CON.y1(1:(end-1-DYN.n_auto),1);                       %solution method vector 
%             obj.mu(1,end+1)         = CON.y1(end,1);                                        %continuation parameter 
%             obj.J(:,:,end+1)        = CON.J1;                                               %Jacobian matrix
%             obj.dy(:,end)           = CON.dy0;                                              %direction vector of the predictor (this one belongs to the previous point)
%             obj.newton_flag(1,1)    = CON.newton_flag;
%             obj.step_width(1,end+1) = CON.step_width;
%             obj.arclength(1,end+1)  = obj.arclength(1,end)+CON.step_width;

        end

        %Abstract Superclass Method: Must be defined

        function IF_arch_bfp_data(obj,CON,DYN,AM)

%             obj.s(:,end+1)          = CON.p_y_bfp(1:(end-1-DYN.n_auto),1);                      %approximation method vector
%             obj.mu(1,end+1)         = CON.p_y_bfp(end,1);                                       %continuation parameter
%             obj.J(:,:,end+1)        = CON.p_J_bfp;                                              %Jacobian matrix
%             obj.dy(:,end)           = NaN.*ones(size(CON.dy0));                                 %direction vector of the predictor (this one belongs to the previous point)
%             obj.step_width(1,end+1) = CON.step_width;
%             obj.arclength(1,end+1)  = CON.p_arclength_bfp;
% 
%             if DYN.n_auto == 0                                                              %Solution is periodic - if this is true: non-autonomous
%                 obj.freq(1,end+1) = DYN.non_auto_freq(CON.p_J_bfp(end,1));
%             elseif DYN.n_auto == 1
%                 obj.freq(1,end+1) = CON.p_J_bfp(end-1,1);
%             end
% 
%             obj.multipliers(:,end+1)    = CON.p_multipliers_bfp;                    %Floquet Multipliers
%             obj.n_unstable(1,end+1)     = 0;                                        %Indiacting number of unstable multipliers
%       
        end


    end

        
    methods(Access = protected) %Can only be called by within class and superclass objects

        %Abstract Superclass Methods: Must be defined with the exact same
        %passed and returned arguments
        
        %Postprocessing
        [s,mu,options]                          = evalsol_solution(obj,DYN,options);            %Function gives back the solution and mu property
        [s_time,mu,time]                        = evalsol_time(obj,DYN,options);                %Function gives back the solution in time domain
        [s_hypertime,mu,hypertimes]             = evalsol_hypertime(obj,DYN,options);           %Function gives back the solution in hypertime domain
        [s_frequency,mu,frequency]              = evalsol_frequency(obj,DYN,options);           %Function gives back the solution in frequency domain

    end

end