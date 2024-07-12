%This is the stability superclass ST. It provides the means for computing
%the (Lyapunov) stability of equilibria, periodic and quasi-periodic solutions
%with different methods
    
    
    classdef Stability < handle
    
        properties
            calc_stability = [];                                            %function handle for the calculation of the stability 
            test_functions = [];                                            %function handle for bifurcation point test functions
            abstol_multiplier = 1e-4;                                       %For bifurcation point iteration: absolute tolerance of the decisive multiplier
            max_iter = 10;                                                  %For bifurcation point iteration: maximum number of iterations
            iterate_bfp = 'on';                                             %Key for iterating the bifurcation point
        end

        properties(GetAccess = public, SetAccess = private)
            curve_container = cell(0);                                      %Container for storing the last n curve points: Needed for iterating a bifurcation point
        end

        properties(Access = private)
            p_container_maxsize = 2;                                           %Number of curve points in the container
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        methods(Access = protected,Sealed)                                  %Protected: Only class and subclass objects can use the method, sealed: method cannot be changed by subclass methods
            [multipliers,vectors,stable] = check_stability_values(obj,multipliers,vectors,stable,stabi_flag); %Method checks the values of the input arguments for inf or NaN values.
            varargout = sort_multipliers(obj,DYN,varargin);                 %General purpose method for sorting multipliers (eigenvalue, Floquet, Lyapunov) and corresponding eigenvectors if possible
            stab_fcn = get_stability_fcn(obj,DYN);                          %When called, this method gives back a vector of the values of a test function (or multipliers) that looses its stability
        end
    
        methods(Static,Sealed)                                                      %Static: Method can be called without creating an object of superclass Stability; Sealed: methods cannot be changed by subclasses
           ST = s_stability_selection(DYN,AM);                                      %Chooses the correct stability subclass and assigns calc_stability according to the solution type and creates it                                                 
           s_ST_gatekeeper(GC,system,opt_stability,opt_sol);                        %Gatekeeper for the Stability class options
        end

        methods
            update_curve_container(obj,DYN,AM,arcl_1,y,multipliers,n_unstable);     %Adds a new curve point and stability information into the curve container and possibly updates the dimension of the contained solution curve points
            update_curve_container_bfp(obj,y_bf,multipliers,n_unstable);            %Updates the curve_container with a curve point from a bifurcation point iteration
            clean_curve_container(obj);                                             %Deletes all elements before the stability change
            y = approx_posc(obj,DYN);                                               %Function approximates the point of stability change based on the curve_container
            [label,msg] = identify_bifurcation(obj);                                %Function give back a string identifying the detected bifurcation
        end

        methods(Abstract, Access = protected) 
            crit_multi(obj,DYN,multipliers);                                                    %Function that needs to be implemented in every subclass: Returns the critical multiplier parts (Real Part, Absolute Value -1, Lyapunov values)
        end

    
    end