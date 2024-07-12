% Class AM_template provides functions which return the residual for shooting
% algorithms. New Method is a AM_template of ApproxMethod
%
%This is a template file

classdef AM_template < ApproxMethod

    properties
        %Define here the properties, which are supplied in the
        %opt_approx_method and opt_init_sol sturctures
        


        %Inherited Properties from Superclass Shoot
        
        %res function_handle
    end
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    methods(Static)                                                        %Static: Method can be called without creating an object of class Shoot;

        s_AM_template_gatekeeper(GC,opt_sol_method,opt_init_sol);          %Gatekeeper method, which is called by the static s_AM_gatekeeper method  to check the inputs at the beginning 

    end
    %%%%%%%%%%%%%%%
    methods
        %% Constructor
        function obj = AM_template(DYN)     
            obj = updateoptions(obj,DYN.opt_approx_method);                %updateoptions method is inherited from ApproxMethod
            obj = updateoptions(obj,DYN.opt_init_sol);                     %updateoptions method is inherited from ApproxMethod
        end
        
        %% Methods
        %Interface methods: This is an abstract method and must be defined
        %Use this interface to pass information between the continuer
        %algorithm and your AM subclass
        function obj = IF_up_res_data(obj,CON)                                     
        
            %Exemplary usage:
%           obj.x0 = CON.yp(1:(end-1));                                            %update the current initial condition. Used for the poincare phase condition.

        end
        
        function obj = getIV(obj,DYN)                                              %Method to generate a method specific initial value from initial condition

            %Exemplary usage:
%           obj.iv = DYN.opt_init.ic;                                               %This is an example for the equilibrium case  

        end

        %% Functions for residuum and computation
        function res = template_residuum_fun(obj,y,DYN)                                   %Method required for shoot_single_nonauto
            
            %Define here (or an external file) the residuum equation. The
            %return value of 

        end

    end
end