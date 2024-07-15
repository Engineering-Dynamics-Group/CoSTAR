% This static function creates an according solution object for the
% specific solution type and solution method
%
%@S:        Solution object
%@DYN:      DynamicalSystem object

function S = s_solution_selection(DYN,AM)


%%%%%%%%%%%%% EQUILIBRIUM %%%%%%%%%%%%%
if(strcmpi(DYN.sol_type,'equilibrium'))
    S = SOL_EQ();


%%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%
elseif(strcmpi(DYN.sol_type,'periodic'))
    
    if(strcmpi(DYN.approx_method,'shooting'))
        S = SOL_PS_SHM(AM);
    end

    if(strcmpi(DYN.approx_method,'mshm'))
        S = SOL_PS_MSHM(AM);
    end
    
    if(strcmpi(DYN.approx_method,'fourier-galerkin'))
        S = SOL_PS_FGM();
    end

    if(strcmpi(DYN.approx_method,'finite-difference'))
        S = SOL_PS_FDM();
    end


%%%%%%%%%%%%% QUASIPERIODIC %%%%%%%%%%%%%
elseif(strcmpi(DYN.sol_type,'quasiperiodic'))

    if(strcmpi(DYN.approx_method,'shooting'))
        S = SOL_QPS_SHM(AM);
    end 

    if(strcmpi(DYN.approx_method,'fourier-galerkin'))
        S = SOL_QPS_FGM();
    end

    if(strcmpi(DYN.approx_method,'finite-difference'))
        S = SOL_QPS_FDM();
    end

end


%Copy the DynamicalSystme ID to the Solution Class id. This helps to match DYN and S object in postprocessing.
S.S_id = DYN.DYN_id;


end