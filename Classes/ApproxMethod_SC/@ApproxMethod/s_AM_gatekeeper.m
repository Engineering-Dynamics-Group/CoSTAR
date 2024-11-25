%Gatekeeper function for the ApproxMethod class. In here, all input parameters are checked, before further processing.
%
%@GC:                   Gatekeeper object
%@system:               User supplied option structure for the system
%@opt_sol:              User supplied option structure for the solution
%@opt_approx_method:    User supplied option structure for the approximation method
%@opt_init:             User supplied option structure for computing the
%initial solution on the curve

function s_AM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init)


%%%%%%%%%%%%% EQUILIBRIUM %%%%%%%%%%%%%
if strcmpi(opt_sol.sol_type,'equilibrium') || strcmpi(opt_sol.sol_type,'eq')
    
        AM_EQ.s_EQ_gatekeeper(GC,system,opt_approx_method,opt_init);


%%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%
elseif strcmpi(opt_sol.sol_type,'periodic') || strcmpi(opt_sol.sol_type,'ps')
    
    if strcmpi(opt_sol.approx_method,'shooting') || strcmpi(opt_sol.approx_method,'shm')
        AM_PS_SHM.s_PS_SHM_gatekeeper(GC,system,opt_approx_method,opt_init);
    end

    if strcmpi(opt_sol.approx_method,'fourier-galerkin') || strcmpi(opt_sol.approx_method,'fgm')
        AM_PS_FGM.s_PS_FGM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init);
    end

    if strcmpi(opt_sol.approx_method,'finite-difference') || strcmpi(opt_sol.approx_method,'fdm')
        AM_PS_FDM.s_PS_FDM_gatekeeper(GC,system,opt_approx_method,opt_init);
    end


%%%%%%%%%%%%% QUASIPERIODIC %%%%%%%%%%%%%
elseif strcmpi(opt_sol.sol_type,'quasiperiodic') || strcmpi(opt_sol.sol_type,'qps')
    
    if strcmpi(opt_sol.approx_method,'fourier-galerkin') || strcmpi(opt_sol.approx_method,'fgm')
        AM_QPS_FGM.s_QPS_FGM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init);
    end

    if strcmpi(opt_sol.approx_method,'shooting') || strcmpi(opt_sol.approx_method,'shm')
        AM_QPS_SHM.s_QPS_SHM_gatekeeper(GC,system,opt_sol,opt_approx_method,opt_init);
    end

    if strcmpi(opt_sol.approx_method,'finite-difference') || strcmpi(opt_sol.approx_method,'fdm')
        AM_QPS_FDM.s_QPS_FDM_gatekeeper(GC,system,opt_approx_method,opt_init);
    end

end


end