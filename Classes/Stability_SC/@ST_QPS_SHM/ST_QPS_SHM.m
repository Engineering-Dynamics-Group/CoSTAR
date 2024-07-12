%This is the stability subclass ST_QPS_SHM. It provides the means for computing
%the (Lyapunov) stability of quasi-periodic solutions via shooting.


classdef ST_QPS_SHM < Stability

    properties

        odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-10);                         % options for the ODE integrators.
        solver string = 'ode45';
        solver_function function_handle
        fsolve_opts = optimoptions('fsolve','Display','none','UseParallel',false,'MaxFunctionEvaluations', 20000,'MaxIter',1e3); %Options for fsolve to solve corrector-equation

        bifurc_label = {'FB','PDB','NSB','BF'};                                                     % Labels for bifurcation, which correspond to the test_functions.
        msg_label = {'This is a fold, pitchfork or transcritical bifurcation.',...                  % The last entry is a safety measure, for the case when a bifurcation
                     'This is a period doubling (aka flip) bifurcation.',...                        % has been detected (due to a change in the number of unstable multipliers),
                     'This is a Neimark-Sacker (aka torus) bifurcation.',...                        % but no test function was active.
                     'A bifurcation has occurred for which a test function does not exist yet.'};

    end


    methods

        function obj = ST_QPS_SHM(DYN)

            obj = updateoptions(obj,DYN.opt_stability);                         % updateoptions method is a general function
                                                                             % Use the solver of the approximation method if given and if not specified in the opt_stability
            if ~isfield(DYN.opt_stability,'solver')
                if isfield(DYN.opt_approx_method,'solver')
                    obj.solver = DYN.opt_approx_method.solver;
                end
            end

            obj = setSolver(obj,obj.solver);                                    % set the solver for the shooting for the monodromy matrix

        end

        [multipliers,vectors,n_unstable,stability_flag] = QPS_SHM_calc_stability(obj,y,J,DYN,AM)     % Method for calculating the Lyapunov exponents
        Output = GSOrthonormalization(obj,MatIn);                               % Gram-Schmidt orthonormalization
        J_korr = jacobi_int(obj,y,AM,Omega,index);                              % Correction Matrix for quasi-periodic stability of shooting method

    end

    methods(Access = protected)
            cm = crit_multi(obj,DYN,multipliers);                                                   %This function returns the absolute value -1  of the given Floquet multipliers. 
    end

    methods(Static)                                                              % Static: Method can be called without creating an object of subclass ST_EQ;
        s_ST_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_stability);                    % Gatekeeper for the Stability subclass options
        help_struct = s_help_opt_stability_QPS_SHM();                                %Help function for the opt_stability struct

    end

end