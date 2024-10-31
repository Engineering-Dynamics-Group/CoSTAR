% Gatekeeper function for the Stability class
% In here, the subsclass gatekeepers are called
%
% @GC:             object of Gatekeeper class
% @system:         user supplied option structure for the system 
% @opt_sol:        user supplied option structure for the solution
% @opt_stability:  user supplied option structure for the stability


function s_ST_gatekeeper(GC,system,opt_sol,opt_stability)

    % Call the subclass gatekeepers

    %%%%%%%%%%%%% EQUILIBRIUM %%%%%%%%%%%%%
    if strcmpi(opt_sol.sol_type,'equilibrium') || strcmpi(opt_sol.sol_type,'eq')

            ST_EQ.s_ST_EQ_gatekeeper(GC,system,opt_sol,opt_stability);


    %%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%
    elseif strcmpi(opt_sol.sol_type,'periodic') || strcmpi(opt_sol.sol_type,'ps')

            ST_PS_SHM.s_ST_PS_SHM_gatekeeper(GC,system,opt_sol,opt_stability);


    % %%%%%%%%%%%%% QUASIPERIODIC %%%%%%%%%%%%%
    elseif strcmpi(opt_sol.sol_type,'quasiperiodic') || strcmpi(opt_sol.sol_type,'qps')

            ST_QPS_SHM.s_ST_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_stability);


    end


end