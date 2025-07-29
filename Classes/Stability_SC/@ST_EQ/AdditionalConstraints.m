%This function provides additional constraints to be checked if a
%bifurcation of an equilibrium solution is detected. This is to check
%whether a fold bifurcation occurs or a transcritcal/pitchfork
%bifurcation occurs instead
%
%@obj:          object of Stability subclass EQ_ST_PS
%@idx0:         Index of bifurcation (1 = FB, 2 = HB)
%@J:            Jacobian at bifurcation point
%
%@idx:          Index for occurring bifurcation (1 = FB, 2 = HB, 3 = TCB/PFB)

function idx = AdditionalConstraints(obj,idx0,J)

n = size(J,1);
dim = n-1;
if(idx0==1)
    s = svd(J(1:end-1,:));                  % Do an SVD of the extended Jacobian [df/ds,df/dmu]
    tol = n*max(abs(s),[],'all')*eps;       % Define tolerance below which value is assumed equal to zero
    s(abs(s) < tol) = 0;                    % Set all values less than tol equal to zero

    rkJext = sum(abs(sign(s)),1);           % Determine the rank of the extended Jacobian

    if(rkJext<dim)
        idx = 3;                            % If rank < dim change index to 3 ('TCB/PFB')
    else
        idx = idx0;                         % If rank == dim leave index at 1 ('FB')
    end
else
    idx = idx0;                             % If index is 2 ('HB') leave index as is
end

end