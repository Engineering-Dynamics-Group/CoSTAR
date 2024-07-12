%This function determines a new higher harmonic and Fourier coefficient based 
%on the residuum in the frequency spectrum. It automatically updates AM.hmatrix
%
%@obj: ApproxMethod AM_PS_FGM object
%@y:   solution approximation continuation vector
%@DYN: DynamicalSystem class object
%
%@yp:  new predicted solution point
%@indisc_flag: 0 if no new higher harmonic could be determined. 1 if everything went alright

function [yp,flag] = IF_increase_discretization(obj,y,DYN)

    hhm0 = obj.hmatrix; %Is needed, if the hmatrix predicted is not unique (see last if - else)
    if (size(obj.hmatrix,2)-1)<obj.n_hh_max %"-1": Constant term is not considered as a higher harmonic
        %Check if the number of evaluation points of the fft needs to be increased in order to be able to predict the necessary higher harmonics
        %This works, since we are only increasing the number of higher harmonics by 1 in every loop.
        if max(obj.hmatrix,[],'all') > obj.n_fft/2
            tmp = nextpow2(2.*max(obj.hmatrix,[],'all')+1);
            warning(append('During increasing the number of higher harmonics in order to meet the error tolerances, the number of FFT evaluations points n_fft was increased from ',num2str(obj.n_fft),' to ',num2str(2^tmp),' in order to comply to the Nyquist-Shannon theorem.'));
            obj.n_fft = 2^tmp;
       end
            
        n0 = obj.n_fft/2+1;
        proj = obj.residuum_projection(y,DYN);

        [~,idx] = max(vecnorm(abs(proj),2,3),[],'all');     %Find the largest element in the projected residuum

        if ~isempty(idx)    %Only move on, if an index has been found

            [r_idx,c_idx] = ind2sub(size(proj,[1,2]),idx);      %Transform the linear index into row and column

            hh(1,1) = r_idx-n0;                                 %Shift the indices
            hh(2,1) = c_idx-n0;

            %max returns the first index of the largest element. Projected residuum is point symmetric. Shift the higher harmonics to the right quadrant.
            if hh(2,1)<0;  hh(1,1)  = - hh(1,1) ; hh(2,1) =abs(hh(2,1)); end
            if hh(2,1) == 0;  hh(1,1)  = abs( hh(1,1) ); end
        
            %Identify the largest element in the projected residuum. Nyquist frequency is considered by obj.n_fft/2-1
            obj.hmatrix = [obj.hmatrix,hh];      %"-1" since, the first element of the projection is the 0th frequency.

        else
            idx = [];
        end
    else
        idx = [];
    end
    
    if isempty(idx)||sum(isnan(obj.hmatrix),'all')||sum(isinf(obj.hmatrix),'all')||~(size(obj.hmatrix,2)==size(unique(obj.hmatrix.','rows').',2)) %Check if an idx was found, there are any inf's or nan's and if the new matrix is unique
        yp = y;
        flag = 0;
        obj.hmatrix = hhm0;
    else
        s0 = y(1:(end-DYN.n_auto-1),1);           %cut off the autonomous frequencies and the curve parameter at the end of the constant vector (stacked Fourier coefficients)
        [s,obj.hmatrix] = obj.sort_guess_FC(DYN,s0);
        yp = [s;y(end-DYN.n_auto:end)];
        flag = 1;
    end

    end