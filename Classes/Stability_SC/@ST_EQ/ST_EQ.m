%This is the stability subclass ST_EQ. It provides the means for computing
%the (Lyapunov) stability of equilibria
    
    
    classdef ST_EQ < Stability
    
        properties
    
            bifurc_label = {'FB','HB','TCB/PFB','BF'};                                                                % Labels for bifurcation, which correspond to the test_functions. 
            msg_label = { 'This is a fold, pitchfork or transcritical bifurcation.',...                     % The last entry is a safety measure, for the case when a bifurcation
                          'This is a Hopf bifurcation.',...                                                 % has been detected (due to a change in the number of unstable multipliers),
                          'A bifurcation has occurred for which a test function does not exist yet.'};      % but no test function was active.
        
        end
    

        methods
            
            function obj = ST_EQ(DYN)
                   obj = updateoptions(obj,DYN.opt_stability);             %updateoptions method is a general function
            end
        
            [multipliers,vectors,stable,max_mult] = EQ_calc_stability(obj,y,J,DYN); %Method for calculating the eigenvalues and eigenvectors of the current solution
            test_fcn = EQ_test_functions(obj,multipliers,DYN);                  %Evaluates the test_functions for determination of the bifurcation type and its position
        
        end
    

        methods(Static)                                                               %Static: Method can be called without creating an object of subclass ST_EQ; 
           s_ST_EQ_gatekeeper(GC,system,opt_sol,opt_stability);                       %Gatekeeper for the Stability subclass options
           help_struct = s_help_opt_stability_EQ();                                   %                               
        end

        
        methods(Access = protected)
            cm = crit_multi(obj,DYN,multipliers);                                       %This function returns the real part of the given eigenvalues. 
            idx = AdditionalConstraints(obj,idx,J);                                     % Checks whether fold-bifurcation is transcritical or pitchfork
        end

    end