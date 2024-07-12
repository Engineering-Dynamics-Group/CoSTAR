% Class ApproxMethod sets the residual which is need for the continuation,
% % according to the chosen Approximation Method (e.g. Shoot). ApproxMethod is a
% superclass

classdef ApproxMethod < handle
    
    properties
        res function_handle                                                 %residual function defined by a solution method (e.g. Shoot)
        iv
        error_control = 'off';                                                  %Switch for using the error_control feature (default value: off)
    end
    %%%%%%%%%%%
    %%%%%%%%%%%
    methods(Static,Sealed)                                                      %Static: Method can be called without creating an object of superclass SolutionType; Sealed: methods cannot be changed by subclasses
        AM = s_method_selection(obj,DYN);                                       %Method for selecting the corresponding Approximation Method (e.g. Shooting) for the Chosen Solution Type (Equilibrium, Periodic, etc.)
        s_AM_gatekeeper(GC,system,opt_sol,opt_sol_method,opt_init);             %Gatekeeper for the ApproxMethod class options
        
    end
    
    methods(Abstract)
        IF_up_res_data(obj);                                                    %The interface update residual data function (IF_up_res_data) is called from the Continuation method to update specific
                                                                                %data for the construction of the residual function. It is extended in the subclasses, if necessary.
        getIV(obj);                                                             %Method to generate a method specific initial value from initial condition
       
        
    end
    
end