%This function determines a new higher harmonic and Fourier coefficient based 
%on the residuum in the frequency spectrum. It automatically updates AM.hmatrix
%
%@obj: ApproxMethod AM_PS_FGM object
%@y1:  solution approximation continuation vector
%@DYN: DynamicalSystem class object
%
%@yp:  new predicted solution point
%@flag: -1 if no new higher harmonic could be determined. 0 if maximal discretization has been reached. 1 if everything went alright.

function [yp,flag] = IF_increase_discretization(obj,y1,DYN)

    hhm0 = obj.hmatrix; %Is needed, if the hmatrix predicted is not unique (see last if - else)

    if (numel(obj.hmatrix)-1)<obj.n_hh_max %"-1": Constant term is not considered as a harmonic

        proj = obj.residuum_projection(y1,DYN);
        [~,idx] = max(vecnorm(abs(2.*proj(:,1:(obj.n_fft/2))),2,1));  %Identify the largest element in the projected residuum. Nyquist frequency is considered by obj.n_fft/2
        obj.hmatrix = [obj.hmatrix,idx-1];      %"-1" since, the first element of the projection is the 0th frequency.
    
        if isempty(idx)||sum(isnan(obj.hmatrix))||sum(isinf(obj.hmatrix))||~(numel(obj.hmatrix)==numel(unique(obj.hmatrix))) %Check if an idx was found, there are any inf's or nan's and if the new matrix is unique
            yp = y1;
            flag = -1;
            obj.hmatrix = hhm0;
        else
            % If we have reached the maximum possible frequency (obj.n_fft/2 is the Nyquist frequency), we need to increase obj.n_fft
            % Theoretically, we could recompute the solution with the updated hmatrix and the old n_fft since we are still below the Nyquist frequency
            % If we did that, however, the estimated error would suddenly become quite small (~ e-11 and below) -> consequence: alternating in- and decrease of discretization and max(hmatrix) does not get higher than n_fft/2-1
            % Recomputing the solution with an increased n_fft leads to a reasonable error again
            if max(obj.hmatrix) == (obj.n_fft/2-1)
                obj.n_fft = 2*obj.n_fft;
                info_text = append('The number of FFT evaluations points n_fft was increased from 2^',num2str(log2(0.5*obj.n_fft)),' to 2^',num2str(log2(obj.n_fft)), ...
                                   ' = ',num2str(obj.n_fft),' to comply with the Nyquist-Shannon criterium and to avoid problems with the error control.');
                write_log(DYN,info_text)                                        % Write info text in log file
                if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                    disp(info_text);                                            % Display information
                end
            end
            s0 = y1(1:(end-DYN.n_auto-1),1);            %cut off the autonomous frequencies and the curve parameter at the end of the constant vector (stacked Fourier coefficients)
            [s,obj.hmatrix] = obj.sort_guess_FC(DYN,s0);
            yp = [s;y1(end-DYN.n_auto:end)];
            flag = 1;
        end
   
    else                    % Maximum number of harmonics has been reached - further increase not possible

        yp = y1;
        flag = 0;
        
    end

end