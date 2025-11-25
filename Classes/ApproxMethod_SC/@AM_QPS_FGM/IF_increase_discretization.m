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
    
    if (size(obj.hmatrix,2)-1)<obj.n_hh_max %"-1": Constant term is not considered as a higher harmonic
            
        n0 = obj.n_fft/2+1;
        proj = obj.residuum_projection(y1,DYN);

        [~,idx] = max(vecnorm(abs(proj),2,3),[],'all');     %Identify the largest element in the projected residuum. Nyquist frequency is considered by obj.n_fft/2-1

        if ~isempty(idx)    %Only move on, if an index has been found

            [r_idx,c_idx] = ind2sub(size(proj,[1,2]),idx);      %Transform the linear index into row and column

            hh(1,1) = r_idx-n0;                                 %Shift the indices
            hh(2,1) = c_idx-n0;

            %max returns the first index of the largest element. Projected residuum is point symmetric. Shift the higher harmonics to the right quadrant.
            if hh(2,1)<0;  hh(1,1)  = - hh(1,1) ; hh(2,1) =abs(hh(2,1)); end
            if hh(2,1) == 0;  hh(1,1)  = abs( hh(1,1) ); end
        
            obj.hmatrix = [obj.hmatrix,hh];      %"-1" since, the first element of the projection is the 0th frequency.

        end

        if isempty(idx)||sum(isnan(obj.hmatrix),'all')||sum(isinf(obj.hmatrix),'all')||~(size(obj.hmatrix,2)==size(unique(obj.hmatrix.','rows').',2)) %Check if an idx was found, there are any inf's or nan's and if the new matrix is unique
            yp = y1;
            flag = -1;
            obj.hmatrix = hhm0;
        else
            % If we have reached the maximum possible frequency (obj.n_fft/2 is the Nyquist frequency), we need to increase obj.n_fft
            % Theoretically, it would be sufficient to increase n_fft/2 if the Nyquist frequency is reached. However, we already increase n_fft when
            % n_fft/2-1 is reached. This is done in the periodic case to avoid some problems and we do it here as well so that the code is identical
            % Moreover: When increasing at n_fft/2, it is very likely that we need another error control iteration since the new error is probably still above the limit due to the increased n_fft value
            if max(abs(obj.hmatrix),[],'all') == (obj.n_fft/2-1)
                obj.n_fft = 2*obj.n_fft;
                info_text = append('The number of FFT evaluations points n_fft was increased from 2^',num2str(log2(0.5*obj.n_fft)),' to 2^',num2str(log2(obj.n_fft)), ...
                                   ' = ',num2str(obj.n_fft),' to comply with the Nyquist-Shannon criterium and to avoid problems with the error control.');
                write_log(DYN,info_text)                                        % Write info text in log file
                if strcmpi(DYN.display,'error-control') || strcmpi(DYN.display,'full')
                    disp(info_text);                                            % Display information
                end
            end
            s0 = y1(1:(end-DYN.n_auto-1),1);           %cut off the autonomous frequencies and the curve parameter at the end of the constant vector (stacked Fourier coefficients)
            [s,obj.hmatrix] = obj.sort_guess_FC(DYN,s0);
            yp = [s;y1(end-DYN.n_auto:end)];
            flag = 1;
        end

    else                    % Maximum number of harmonics has been reached - further increase not possible

        yp = y1;
        flag = 0;
        
    end
    
end