% Class Equilibrium provides the residual function for a equilibrium
% continuation. Equilibrium is a subclass of ApproxMethod

classdef AM_EQ < ApproxMethod

    properties
        % Only inherited properties
        
    end

    %%%%%%%%%%%%%%%

    methods(Static)                                                     % Static: Method can be called without creating an object of class AM_EQ

        s_EQ_gatekeeper(GC,system,opt_sol_method,opt_init);             % Gatekeeper method, which is called by the static ST_gatekeeper method 
        help_text = s_help_opt_approx_method_EQ();                      % Help file for the approx_method option structure 
        help_text = s_help_opt_init_EQ();                               % Help file for the approx_method option structure 
    
    end

    %%%%%%%%%%%%%%%
        
    methods 

        % Constructor
        function obj = AM_EQ(DYN)
            
            obj = obj.getIV(DYN);
            
        end
        
        % Functions for equilibrium algorithm
        function [res,J_res] = residual_function(obj,y,DYN)             % Set up the residual function
                
                %For some reason it is way faster to preallocate the variables first and then use them... 
                % it could have something to do with Matlab intern code optimization
                Fcn = DYN.rhs;                                          % Get residual function
                x = y(1:(end-1));                                       % Get solution vector (without continuation parameter)
                mu = y(end);                                            % Get continuation parameter
                
                %Evaluate the active parameter 
                param = DYN.param;
                param{DYN.act_param} = mu;

                res = Fcn(x,param);                                     % Define residual function  

                % Calculate the Jacobian
                dFcn_dz = DYN.jacobian(x,param);                        % df/dz
                h = eps^(1/3)*(1+abs(mu));                              % Finite difference step width
                param_mu_plus = param;  param_mu_minus = param;
                param_mu_plus{DYN.act_param} = mu + h;                  % param array with "+ h" perturbed mu-value
                param_mu_minus{DYN.act_param} = mu - h;                 % param array with "- h" perturbed mu-value
                dFcn_dmu = (Fcn(x,param_mu_plus) - Fcn(x,param_mu_minus)) ./ (2*h);     % df/dmu
                J_res = [dFcn_dz, dFcn_dmu];                            % Jacobian

        end

        % Function wrapper that sets the complete residuum function for the initial solution
        function [F,J] = res_fun_init(obj,y,y0)

            [res,J_res] = obj.res(y);
            F = [res; y(end)-y0(end)];
            J = [J_res; zeros(1,length(y)-1), 1];

        end

        % Function wrapper that sets the complete residuum function
        function [F,J] = res_fun(obj,y,CON)

            [res,J_res] = obj.res(y);
            F = [res; CON.sub_con(y,CON)];
            J = [J_res; CON.d_sub_con(y,CON)];

        end

        % Interface methods
        function obj = IF_up_res_data(obj,var,DYN); end                 % Nothing needs to be done here
        
        obj = getIV(obj,DYN);                                           % Compute initial value

        % updateoptions (not used here)

    end

end