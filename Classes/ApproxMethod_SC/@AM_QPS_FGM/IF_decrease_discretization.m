%This function decreases the descretization of the Fourier series.
%For smooth systems, the highest harmonics has always the lowest impact on the error. However this might not be the case for non-smooth systems, which
%is why we test, which harmonic has the lowest impact
%
%@obj: ApproxMethod AM_QPS_FGM object
%@y1:  solution approximation continuation vector
%@DYN: DynamicalSystem class object
%
%@yp:  new predicted solution point
%@flag: 1 if everything went alright. 3 if something went wrong. -1 if minimum number of harmonics has already been reached

function [yp,flag] = IF_decrease_discretization(obj,y1,DYN)
    
    hhm0 = obj.hmatrix;                       %store the current value for later use
    
    %Storing the hmatrix to an additional property: 
    %This is needed for the IF_update_sol_dim method to have a reference point in order to determine which higher harmonic was eliminated
    obj.p_hmatrix_old = obj.hmatrix;    

    %reshape Fourier coefficients
    tmp = reshape(y1(1:end-DYN.n_auto-1,1),DYN.dim,[]);
    c0 = tmp(:,1);
    cmatrix = tmp(:,2:0.5*(numel(tmp)/DYN.dim+1));          %Gatekeeper has ensured that this computation of dimension is always possible
    smatrix = tmp(:,0.5*(numel(tmp)/DYN.dim+1)+1:end);
    

    if size(obj.hmatrix,2)>3     %lowest number for higher harmonics

        err = zeros(1,size(hhm0,2)-3);
        for k =3:(size(hhm0,2)-1)    %Delete iteratively one frequency to find the frequency with the lowest impact on the error, which is then deleted.
            idx = [1:(k-1),(k+1):size(hhm0,2)-1];
            obj.hmatrix = [zeros(2,1),hhm0(:,idx+1)];
            ytmp = [c0; reshape(cmatrix(:,idx),[],1); reshape(smatrix(:,idx),[],1);y1(end-DYN.n_auto:end)];
            err(1,k-2) = obj.IF_estimate_error(ytmp,DYN);
        end
    
        obj.hmatrix = hhm0;         %Set hmatrix to the initial value
        [~,idx_min] = min(err);     %idx+2 is the harmonic with the smallest impact on the error, which is consequently delted
        idx_min = idx_min + 2;

        if isempty(idx_min)     %Check if an idx was found, there are any inf's or nan's and if the new matrix is unique
            yp = y1;
            flag = 3;
        else                    % Reduce the number of harmonics by omitting the harmonic with the smallest impact on the error
            idx = [1:(idx_min-1),(idx_min+1):size(hhm0,2)-1];
            obj.hmatrix = [zeros(2,1),hhm0(:,idx+1)];
            yp = [c0; reshape(cmatrix(:,idx),[],1); reshape(smatrix(:,idx),[],1);y1(end-DYN.n_auto:end)];
            flag = 1;
            % If max(hmatrix) falls below n_fft/4-1, we can halve n_fft
            % Reason: If n_fft is set to n_fft/2, the highest possible harmonic according to Nyquist-Shannon is (n_fft/2)/2-1 = n_fft/4-1. Additionally, we have to consider that n_fft ...
            %         is automatically doubled in IF_increase_discretization if the highest possible harmonic is reached. Ergo, max(hmatrix) < n_fft/4-1 to reasonable halve n_fft
            % Apart from that, we also have to consider that the minimum value of n_fft is either the default value or the value set by the user
            if isfield(DYN.opt_approx_method,'n_fft');  n_fft_min = DYN.opt_approx_method.n_fft;        % Value set by user
            else;                                       n_fft_min = 2^6;                                % Default value
            end
            if (max(obj.hmatrix,[],'all') <= (obj.n_fft/4-2)) && (obj.n_fft > n_fft_min)    % Halve n_fft in accordance with Nyquist-Shannon and if current n_fft is larger than minimum n_fft
                obj.n_fft = 0.5*obj.n_fft;
                info_text = append('The number of FFT evaluations points n_fft was decreased from 2^',num2str(log2(2*obj.n_fft)),' to 2^',...
                                    num2str(log2(obj.n_fft)),' = ',num2str(obj.n_fft),' in accordance with the Nyquist-Shannon criterium to reduce computing effort.');
                write_log(DYN,info_text)                                        % Write info text in log file
                if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                    disp(info_text);                                            % Display information
                end
            end
        end

    else                        % Minimum number of harmonics has been reached - further reduce not possible

        yp = y1;
        flag = 2;

    end
   
end