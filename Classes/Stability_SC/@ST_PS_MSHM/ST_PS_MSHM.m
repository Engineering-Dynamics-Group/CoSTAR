    %This is the stability subclass ST_PS_SHM. It provides the means for computing
    %the (Lyapunov) stability of periodic solutions via shooting.
    
    
    classdef ST_PS_MSHM < Stability
    
        properties
    
            odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);                    %options for the ODE integrators.
            solver string = 'ode45';
            n_shoot = 5;
            solver_function function_handle
            fsolve_opts = optimoptions('fsolve','Display','none','UseParallel',false,'MaxFunctionEvaluations', 20000,'MaxIter',1e3); %Options for fsolve to solve corrector-equation

            bifurc_label = {'FB','PDB','NSB','BF'};                                             % Labels for bifurcation, which correspond to the test_functions.
            msg_label = {'This is a fold, pitchfork or transcritical bifurcation.',...          % The last entry is a safety measure, for the case when a bifurcation
                         'This is a period doubling (aka flip) bifurcation.',...                         % has been detected (due to a change in the number of unstable multipliers),
                         'This is a Neimark-Sacker (aka torus) bifurcation.',...                         % but no test function was active.
                         'A bifurcation has occurred for which a test function does not exist yet.'};

        end
    
        methods
    
            function obj = ST_PS_MSHM(DYN)
    
                obj = updateoptions(obj,DYN.opt_stability);             %updateoptions method is a general function
    
                %Use the solver of the approximation method if given and if not specified in the opt_stability
                if ~isfield(DYN.opt_stability,'solver')
                    if isfield(DYN.opt_approx_method,'solver')
                        obj.solver = DYN.opt_approx_method.solver;
                    end
                end
                obj = setSolver(obj,obj.solver);                        %set the solver for the shooting for the monodromy matrix
    
            end
    
            [multipliers,vectors,stable,max_mult] = PS_SHM_calc_stability_non_auto(obj,y,J,DYN,AM);        %Function for caluclating Floquet multipliers for a non-autonomous periodic solution for different approximation methods
            [multipliers,vectors,stable,max_mult] = PS_SHM_calc_stability_auto(obj,y,J,DYN,AM);            %Function for caluclating Floquet multipliers for an autonomous periodic solution for different approximation methods
     

            res = SHM_single_fun(obj,x,mu,DYN);                                                   %ODE integration function for use in a shooting method for non-autonomous systems: Here, for usage in computing Floquet multipliers 
            res = SHM_single_auto_fun(obj,x,mu,x0,DYN);                                            %ODE integration function for use in a shooting method for autonomous systems: Here, for usage in computing Floquet multipliers
            
            res = MSHM_fun(obj,x,mu,DYN);
            res = MSHM_auto_fun(obj,x,mu,DYN);

        end
    
        methods(Access = protected)
            cm = crit_multi(obj,DYN,multipliers);                                                   %This function returns the absolute value -1  of the given Floquet multipliers. 
        end

        methods(Static)                                                                 %Static: Method can be called without creating an object of subclass ST_PS_SHM;
            s_ST_PS_SHM_gatekeeper(GC,system,opt_sol,opt_stability);                        %Gatekeeper for the Stability subclass options
            help_struct = s_help_opt_stability_PS_SHM();                                %Help function for the opt_stability struct
            tfcn = PS_test_functions(multipliers,DYN);                                  %Test function for periodic solutions.
        end
    
    end