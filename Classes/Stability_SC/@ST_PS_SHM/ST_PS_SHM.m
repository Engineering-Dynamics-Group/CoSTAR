% This is the stability subclass ST_PS_SHM. It provides the means for computing
% the (Lyapunov) stability of periodic solutions via shooting.


classdef ST_PS_SHM < Stability

    properties

        odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);                 % options for the ODE integrator
        solver string = 'ode45';
        n_shoot = 5;
        solver_function function_handle
        fsolve_opts = optimoptions('fsolve','Display','none','MaxIter',50,'UseParallel',false,'SpecifyObjectiveGradient',true);    % Options for fsolve to solve corrector-equation

        bifurc_label = {'FB','PDB','NSB','BF'};                                             % Labels for bifurcation, which correspond to the test_functions.
        msg_label = {'This is a fold, pitchfork or transcritical bifurcation.',...          % The last entry is a safety measure, for the case when a bifurcation
            'This is a period doubling (aka flip) bifurcation.',...                         % has been detected (due to a change in the number of unstable multipliers),
            'This is a Neimark-Sacker (aka torus) bifurcation.',...                         % but no test function was active.
            'A bifurcation has occurred for which a test function does not exist yet.'};

    end

    %%%%%%%%%%%%%%%

    methods

        % Constructor
        function obj = ST_PS_SHM(DYN)

            % Reset the default value for the solver if shooting is the approximation method and solver was not defined in opt_stability
            if ~isfield(DYN.opt_stability,'solver')
                if isfield(DYN.opt_approx_method,'solver')
                    obj.solver = DYN.opt_approx_method.solver;
                end
            end
            obj = setSolver(obj,obj.solver);                            % Set the solver for the shooting for the monodromy matrix
            
            obj = updateoptions(obj,DYN.opt_stability);                 % updateoptions method is a general function

        end

        [res,J_res] = PS_SHM_ST_residuum(obj,x,x0,mu,DYN);              % Residuum function
        [multipliers,vectors,stable,max_mult] = PS_SHM_calc_stability(obj,y,J,DYN,AM);      % Function for calculating Floquet multipliers for different approximation methods
        
    end

    %%%%%%%%%%%%%%%

    methods(Access = protected)

        cm = crit_multi(obj,DYN,multipliers);                           % This function returns the absolute value -1  of the given Floquet multipliers
        idx = AdditionalConstraints(obj,idx,J);                         % Checks additional crtiteria for bifurcations

    end

    %%%%%%%%%%%%%%%

    
    methods(Static)                                                     % Static: Method can be called without creating an object of subclass ST_PS_SHM
        
        s_ST_PS_SHM_gatekeeper(GC,system,opt_sol,opt_stability);        % Gatekeeper for the Stability subclass options
        help_struct = s_help_opt_stability_PS_SHM();                    % Help function for the opt_stability struct
        tfcn = PS_test_functions(multipliers,DYN);                      % Test function for periodic solutions

    end


end