% Class Equilibrium provides the residual function for a equilibrium
% continuation. Equilibrium is a subclass of ApproxMethod
classdef AM_Equilibrium < ApproxMethod

    properties
        %Inherited properties
        
    end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% 
    methods(Static)                                                        %Static: Method can be called without creating an object of class Equilibrium;

        s_EQ_gatekeeper(GC,system,opt_sol_method,opt_init);                %Gatekeeper method, which is called by the static ST_gatekeeper method 
        help_text = s_help_opt_approx_method_EQ();                         %Help file for the approx_method option structure 
        help_text = s_help_opt_init_EQ();                                  %Help file for the approx_method option structure 
    end
    %%%%%%%%%%%%%%%
    
    
    methods 
        function obj = AM_Equilibrium(DYN)
            
            obj = obj.getIV(DYN);
            
        end
        
        %% Functions for equilibrium algorithm
        function res = residual_function(obj,y,DYN)                        %Set up the residual function
                
                %For some reason it is way faster to preallocate the
                %variables first and then use them... it could have s.th.
                %to do with Matlab intern code optimization
                Fcn = DYN.rhs;                                             %Get residual function
                x = y(1:(end-1));                                          %Get solution vector (without continuation parameter)
                mu = y(end);                                               %Get continuation parameter
                
                %Evaluate the active parameter 
                param = DYN.param;
                param{DYN.act_param} = mu;

                res = Fcn(x,param);                                        %Define residual function  
        end

        %Abstract methods:
        function f = IF_up_res_data(obj,CON); end                          %Nothing needs to be done here
        
        obj = getIV(obj,DYN);                                               %Compute initial value
        %updateoptions (not used here)

    end

end