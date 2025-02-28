%This is the stability subclass ST_QPS_SHM. It provides the means for computing
%the (Lyapunov) stability of quasi-periodic solutions via shooting.


classdef ST_QPS_SHM < Stability

    properties

        n_char_st = 100;                                                % Number of characteristics used to interpolate the mapping matrices
        n_map     = 2e4;                                                % Number of mappings to compute the Ljapunov exponent of m-th order
        solver string = 'ode45';
        solver_function function_handle
        fsolve_opts = optimoptions('fsolve','Display','none','UseParallel',false,'MaxFunctionEvaluations',20000,'MaxIter',50,'SpecifyObjectiveGradient',true); %Options for fsolve to solve corrector-equation

        bifurc_label = {'FB','PDB','NSB','BF'};                                                     % Labels for bifurcation, which correspond to the test_functions.
        msg_label = {'This is a fold, pitchfork or transcritical bifurcation.',...                  % The last entry is a safety measure, for the case when a bifurcation
                     'This is a period doubling (aka flip) bifurcation.',...                        % has been detected (due to a change in the number of unstable multipliers),
                     'This is a Neimark-Sacker (aka torus) bifurcation.',...                        % but no test function was active.
                     'A bifurcation has occurred for which a test function does not exist yet.'};

        QPS_SHM_ST_residuum function_handle                             % Residuum function for reshooting the solution
        AM_ST                                                           % Object of AM_QPS_SHM used to reshoot a solution. 

    end


    methods

        % Constructor
        function obj = ST_QPS_SHM(DYN,AM)

            % Reset the default value for the solver if shooting is the approximation method and solver was not defined in opt_stability
            if ~isfield(DYN.opt_stability,'solver')                     
                if isfield(DYN.opt_approx_method,'solver')              % Use the solver from shooting (approximation method)
                    obj.solver = DYN.opt_approx_method.solver;          % If solver was not defined in opt_approx_method: default value from here is used
                end                                                     
            end
            obj = setSolver(obj,obj.solver);                            % Set the solver for calculating the characteristics

            % Reset the default value for the number of characteristics if shooting is the approximation method and n_char_st was not defined in opt_stability
            if ~isfield(DYN.opt_stability,'n_char_st')                  
                if isfield(DYN.opt_approx_method,'n_char')
                    obj.n_char_st = DYN.opt_approx_method.n_char;       % Use the solver of shooting approximation method
                end
            end

            % Update the properties from the options
            obj = updateoptions(obj,DYN.opt_stability);                 % updateoptions is a general function

            % Use the correct residuum function
            if DYN.n_auto==0
                obj.QPS_SHM_ST_residuum = @(x,mu) qp_SHM_ST_non_auto_fun(obj,x,mu,DYN);             % Residual for non-autonomous case
            elseif DYN.n_auto==1
                obj.QPS_SHM_ST_residuum = @(x,mu) qp_SHM_ST_mixed_fun(obj,x,mu,DYN);                % Residual for mixed case
            elseif DYN.n_auto==2
                obj.QPS_SHM_ST_residuum = @(x,mu) qp_SHM_ST_auto_fun(obj,x,mu,DYN);                 % Residual for full autonomous case
            end

            % Create an object of AM_QPS_SHM if the solution needs to be reshooted
            if ~(strcmpi(DYN.approx_method,'shooting') && (AM.n_char == obj.n_char_st))
                options.system = DYN.system;                                                        % Copy system options
                options.opt_sol = DYN.opt_sol;                                                      % Copy solution options
                options.opt_init = struct();                                                        % Can be empty - obj.iv is filled later on in QPS_SHM_calc_stability by IF_up_res_data
                options.opt_approx_method = struct('solver',obj.solver,'n_char',obj.n_char_st);     % Set solver and n_char using values from here
                options.opt_cont = DYN.opt_cont;                                                    % Copy continuation options
                options.opt_stability = DYN.opt_stability;                                          % Copy stability options
                DYN_ST = DynamicalSystem(options);                                                  % Create new DynamicalSystem object (only needed to create AM_QPS_SHM object)
                obj.AM_ST = AM_QPS_SHM(DYN_ST);                                                     % Create object of AM_QPS_SHM and store the handle in obj.AM_ST
            end

        end

        % Other methods
        [f,J] = qp_SHM_ST_non_auto_fun(obj,x,mu,DYN);                   % Residual for non-autonomous case for quasi-periodic shooting algorithm
        [f,J] = qp_SHM_ST_mixed_fun(obj,x,mu,DYN);                      % Residual for mixed case for quasi-periodic shooting algorithm
        [f,J] = qp_SHM_ST_auto_fun(obj,x,mu,DYN);                       % Residual for full-autonomous case for quasi-periodic shooting algorithm

        [multipliers,vectors,n_unstable,stability_flag] = QPS_SHM_calc_stability(obj,y,J,DYN,AM)     % Method for calculating the Lyapunov exponents
        Output = GSOrthonormalization(obj,MatIn);                       % Gram-Schmidt orthonormalization
        J_korr = jacobi_int(obj,y,AM,Omega,index);                      % Correction Matrix for quasi-periodic stability of shooting method

    end


    methods(Access = protected)
            cm = crit_multi(obj,DYN,multipliers);                       % This function returns the absolute value -1  of the given Floquet multipliers. 
    end


    methods(Static)                                                     % Static: Method can be called without creating an object of subclass ST_EQ;
        s_ST_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_stability);       % Gatekeeper for the Stability subclass options
        help_struct = s_help_opt_stability_QPS_SHM();                   % Help function for the opt_stability struct

    end


end