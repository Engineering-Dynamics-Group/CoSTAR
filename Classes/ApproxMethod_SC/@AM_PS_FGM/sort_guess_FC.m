%  If there are less Fourier coefficients than higher harmonics in obj.hmatrix,
% additional Fourier coefficients are guessed. There are no real limitations on the initial ordering of obj.hmatrix.
% the only requirement is, that the Fourier coefficients correspond to the positions
% in the hmatrix, which also might contain more frequencies (which must not necessarily be "higher").
%
% Example: C0 = zeros(2,1), Cmatrix = Smatrix = [1,2;3,4], hmatrix = [0,1,3,7,2,5,11,9]. 
%
%@obj:          ApproxMethod class object of AM_PS_FGM
%@DYN:          DynamicalSystem class object
%@FC0:          Initial Fourier Coefficient vector
%
%@s:            update solution vector in approximation space (without auto-frequencies or continuation parameter)
%@hmatrix:      update matrix of higher harmonics

function [s,hmatrix] = sort_guess_FC(obj,DYN,FC0)
        
    %Either one of the two is present (ensured by the gatekeeper)
    
    tmp = reshape(FC0,DYN.dim,[]);
    c0 = tmp(:,1);
    cmatrix = tmp(:,2:0.5*(numel(tmp)/DYN.dim+1));          %Gatekeeper has ensured that this computation of dimension is always possible
    smatrix = tmp(:,0.5*(numel(tmp)/DYN.dim+1)+1:end);
    
    %Sort the hmatrix and the FCs accordingly (also account for case when there are more higher harmonics than Fourier coefficients)
    idx_tmp = min([size(cmatrix,2),numel(obj.hmatrix)-1]);                 %-1 due to constant higher harmonic
    hhm = obj.hmatrix;
    
    hhm(hhm==0) = [];                                                          %delete the zero frequency for now. Needed for sorting so that the indices of sorted hmatrix can be applied to cmatrix and smatrix
    [hh_tmp,idx_tmp] = sort(hhm(1:idx_tmp));                                  
    
    cmatrix = cmatrix(:,idx_tmp);
    smatrix = smatrix(:,idx_tmp);
    
    
    if size(cmatrix,2)==(obj.p_n_hh()-1)
        s = [c0(:);cmatrix(:);smatrix(:)];         %For the case that obj.fc0 is given for all higher harmonics in hmatrix... this is a stupid way to do it. But for any other case it is easier to implement.
        hmatrix = [0,hh_tmp];
    else
        % We are guessing the additionally needed initial conditions based on the knowledge, that Fourier coefficients drop exponentially:
        % C_n = K/(n^m). K is always C_1.
        hh_guess = sort(setdiff(hhm,hh_tmp));   %Get the higher harmonics which we want to guess. 
    
        numb_guess = numel(hh_guess);  %Number of needed guesses
        if size(cmatrix,2)>1
            %         m can be determined by the following formula:
            m_c = log2(abs(cmatrix(:,1))./abs(cmatrix(:,end)))./log2(size(cmatrix,2));
            m_s = log2(abs(smatrix(:,1))./abs(smatrix(:,end)))./log2(size(smatrix,2));
    
            %Catch any inf's or nan's
                [ii,jj]=find(isnan(m_c)| isinf(m_c));
                m_c(ii,jj) = 1; 
                [ii,jj]=find(isnan(m_s)| isinf(m_s));
                m_s(ii,jj) = 1; 
      
        else 
            %If we only have the coefficients for the first frequency, 
            % the coefficients cannot be determined by that data alone and we guess them to 
            %be 1 
            m_c = ones(DYN.dim,numb_guess);
            m_s = ones(DYN.dim,numb_guess);
    
        end
    
        %Guess the additional initial conditions.
        tmp_cmatrix =  [cmatrix,repmat(cmatrix(:,1),[1,numb_guess])./hh_guess.^m_c];        % abs() around cmatrix removed
        tmp_smatrix =  [smatrix,repmat(smatrix(:,1),[1,numb_guess])./hh_guess.^m_s];        % abs() around smatrix removed
    
        %Sort the complete h-matrix again. Why? Maybe ferquencies needed to be guessed which are not higher than the hh_tmp frequencies but in between (unlikely case... but now accounted for).
        [hhm,idx_tmp] = sort([hh_tmp,hh_guess]);  
        hmatrix = [0,hhm];
        
        cmatrix = tmp_cmatrix(:,idx_tmp);
        smatrix = tmp_smatrix(:,idx_tmp);
        
        s = [c0(:);cmatrix(:); smatrix(:)];
    
    end




end