    %This function decreases the descretization of the Fourier series.
    %For smooth systems, the highest harmonics has always the lowest impact on the error. However this might not be the case for non-smooth systems, which
    %is why we test, which harmonic has the lowest impact
    %
    %@obj: ApproxMethod AM_QPS_FGM object
    %@y:   solution approximation continuation vector
    %@DYN: DynamicalSystem class object
    %
    %@yp:  new predicted solution point
    %@flag: 1 if everything went alright. 0 if something went wrong. -1 if minimum number of harmonics has already been reached
    
    function [yp,flag] = IF_decrease_discretization(obj,y,DYN)
    
    hhm0 = obj.hmatrix;                       %store the current value for later use
    
    %Storing the hmatrix to an additional property: This is needed for the IF_update_sol_dim method to have a reference point, which higher harmonics were
    %eliminated. 
    obj.p_hmatrix_old = obj.hmatrix;    

    %reshape Fourier coefficients
    tmp = reshape(y(1:end-DYN.n_auto-1,1),DYN.dim,[]);
    c0 = tmp(:,1);
    cmatrix = tmp(:,2:0.5*(numel(tmp)/DYN.dim+1));          %Gatekeeper has ensured that this computation of dimension is always possible
    smatrix = tmp(:,0.5*(numel(tmp)/DYN.dim+1)+1:end);
    
    if size(obj.hmatrix,2)>3     %lowest number for higher harmonics
        err = zeros(1,size(hhm0,2)-3);
        for k =3:(size(hhm0,2)-1)    %Delete iteratively one frequency to find the frequency with the lowest impact on the error, which is then deleted.
            idx = [1:(k-1),(k+1):size(hhm0,2)-1];
            obj.hmatrix = [zeros(2,1),hhm0(:,idx+1)];
            ytmp = [c0; reshape(cmatrix(:,idx),[],1); reshape(smatrix(:,idx),[],1);y(end-DYN.n_auto:end)];
            err(1,k-2) = obj.IF_estimate_error(ytmp,DYN);
        end
    
        obj.hmatrix = hhm0;         %Set hmatrix to the initial value
        [~,idx_min] = min(err);     %idx+1 is the harmonic with the smallest impact on the error, which is consequently delted
        idx_min = idx_min + 2;

        if isempty(idx_min)     %Check if an idx was found, there are any inf's or nan's and if the new matrix is unique
            yp = y;
            flag = 0;
        else                    % Reduce the number of harmonics by omitting the harmonic with the smallest impact on the error
            idx = [1:(idx_min-1),(idx_min+1):size(hhm0,2)-1];
            obj.hmatrix = [zeros(2,1),hhm0(:,idx+1)];
            yp = [c0; reshape(cmatrix(:,idx),[],1); reshape(smatrix(:,idx),[],1);y(end-DYN.n_auto:end)];
            flag = 1;
        end

    else                        % Minimum number of harmonics has been reached - further reduce not possible
        yp = y;
        flag = -1;              % Set special exit flag

    end
   
    end