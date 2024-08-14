% Function defines the residual whoose roots define the curve which is to
% be continued according to the chosen approximation method (e.g. SHM)

function AM = s_method_selection(DYN)


%%%%%%%%%%%%% EQUILIBRIUM %%%%%%%%%%%%%
if(strcmpi(DYN.sol_type,'equilibrium'))
        AM = AM_Equilibrium(DYN);
        AM.res = @(y) AM.residual_function(y,DYN);                          %This call could also be integrated into the constructor of the class - not done for uniformity of this script


%%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%
elseif(strcmpi(DYN.sol_type,'periodic'))
    
    if(strcmpi(DYN.approx_method,'shooting'))
        AM = AM_PS_MSHM(DYN);
        if(DYN.n_auto==0)
                AM.res = @(y)AM.MSHM_fun(y,DYN);                    %Residual for non-autonomous single shooting
        elseif(DYN.n_auto==1)
                AM.res = @(y)AM.MSHM_auto_fun(y,DYN);               %Residual for autonomous single shooting
        end
    end

    if(strcmpi(DYN.approx_method,'fourier-galerkin'))
        AM = AM_PS_FGM(DYN);
        AM.res = @(y)AM.PS_FGM_residuum(y,DYN);  
    end

    if(strcmpi(DYN.approx_method,'finite-difference'))
        AM = AM_PS_FDM(DYN);
        AM.res = @(y) AM.PS_FDM_residuum(y,DYN);                            %Residual function
    end


%%%%%%%%%%%%% QUASIPERIODIC %%%%%%%%%%%%%
elseif(strcmpi(DYN.sol_type,'quasiperiodic'))

    if(strcmpi(DYN.approx_method,'fourier-galerkin'))
        AM = AM_QPS_FGM(DYN);
        AM.res = @(y)AM.QPS_FGM_residuum(y,DYN);  
    end
    
    if(strcmpi(DYN.approx_method,'shooting'))
        AM = AM_QPS_SHM(DYN);
        if(DYN.n_auto==0)                                                   %Residual for non-autonomous case
            AM.res = @(y,DYN)AM.qp_SHM_non_auto_fun(y,DYN);
        elseif(DYN.n_auto==1)                                               %Residual for mixed case
            AM.res = @(y,DYN)AM.qp_SHM_mixed_fun(y,DYN);
        elseif(DYN.n_auto==2)                                               %Residual for full autonomous case
            AM.res = @(y,DYN)AM.qp_SHM_auto_fun(y,DYN); 
        end
    end

    if(strcmpi(DYN.approx_method,'finite-difference'))
        AM = AM_QPS_FDM(DYN);
        AM.res = @(y) AM.QPS_FDM_residuum(y,DYN);                           %Residual function
    end

end


end