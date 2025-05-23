% Class Solution creates an object in which all relevant information for
% the continuation and the solution of the continuation is stored.
% This is a superclass

classdef Solution < handle

    %%%%%%%%%%%%%% Properties %%%%%%%%%%%%%%%

    properties

        y0                      %initial curve point vector
        s                       %storing method solution vector. Undecided if 2D matrix or cell.
        J                       %storing Jacobian matrices. Undecided if 3D matrix or cell.
        mu                      %storing continuation parameter. array
        dy                      %storing the direction vector of the predictor (tangent, secant, ...)
        error                   %solution type specific error
        newton_flag             %storing newton-solver flags. 2D matrix or cell
        stability_flag          %storing stability computation exitflag
        stopping_flag           %information on the cause of ending the computation:

        step_width              %Step width for each solution point
        arclength               %Arclength of Solution curve

        n_unstable              %Indicates the number of unstable multipliers
        multipliers             %Multipliers to determine stability (e.g. Eigenvalues, Floquet-Multipliers or Lyapunov-Exponents)
        max_multiplier          %Maximal multiplier indicating the stability loss. Absolute value of Floquet multipliers are centered around zero.
        vectors                 %(Eigen-)vectors corresponding to multipliers

        bifurcation = array2table(zeros(0,3)); %Empty Table for storing bifurcation information into it (gets its property names in the constructor)
        warnings = cell(0,0);   %Cell array to store warnings

        S_id                    %Solution class identification string 
                                %(equal to the ID of DynamicalSystem (DYN_id) to which the Solution belongs, used in postprocessing to ensure DYN and S match)
    end


    properties(Hidden = true)

        plot_color              %Struct for different standard plotting colors

    end


    %%%%%%%%%%%%%% Methods %%%%%%%%%%%%%%%

    methods

        %Constructor
        function obj = Solution()
            obj.set_plot_colors();                                                         %setting the plot color for plotting in contplot and solplot
            obj.bifurcation.Properties.VariableNames = {'bifurcation','index','info'};     %Set the default property names for the table
        end

        %Postprocessing functions
        [s,mu,varargout]    = solget(obj,DYN,options);                  %Extracts a solution in either solution, time, frequency or torus domain.
        varargout           = solplot(obj,DYN,options);                 %Plots single solution for specific datapoint(s)
        varargout           = contplot(obj,DYN,options);                %Plots continued curve

    end


    methods(Static,Sealed)

        S = s_solution_selection(obj,DYN);      %Creates the according solution object for the solution type and method chosen.
        help_struct = s_help_contplot();        %costarhelp function for the contplot options
        help_struct = s_help_solplot();         %costarhelp function for the solplot options
        help_struct = s_help_solget();          %costarhelp function for the solget options

    end


    methods(Access = protected, Sealed) %Only accessed by methods of the class Solution

        check_fcn_handle(obj,DYN,GC,fcn_handle_options,name_handle,d_in,d_out)                      %Small function checking, if the function_handle is correct

        options = solget_gatekeeper(obj,DYN,options);                                               %Checks the options structure entries for solget
        options = solplot_gatekeeper(obj,DYN,options);                                              %Checks the options structure entries for solplot
        options = contplot_gatekeeper(obj,DYN,options);                                             %Checks the options structure entries for contplot

        options = solget_up_index(obj,DYN,options);                                                 %Updates options.index for the evaluation of the indices
        set_plot_colors(obj);                                                                       %sets the plot_color property of the Solution class object
        rgb = custom_color(obj,options);                                                            %Gives either back a User defined or random color rgb code

    end


    methods(Abstract)                           %Abstract methods must be implemented in the subclasses. Additional in- and output can be added

        IF_arch_init_data(obj);                 %Interface method for archiving the data of the initial solution
        IF_arch_data(obj);                      %Interface method for archiving the data of the continuation
        IF_arch_bfp_data(obj,CON,DYN,AM);       %Interface method for archiving the data of a found bifurcation point into the curve

    end


    methods(Abstract, Access = protected)       %Only class and subclass can access these methods

        evalsol_time(obj,DYN,options);          %Function gives back the solution in time domain
        evalsol_hypertime(obj,DYN,options);     %Function gives back the solution in hypertime space
        evalsol_frequency(obj,DYN,options);     %Function gives back the solution in frequency domain


    end


end