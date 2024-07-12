%   This is a method of the Stability subclass ST_QPS_SHM
%   It calculates the Lyapunov exponents and indicates if the solution is stable or unstable
%
%@obj:  ST_QPS_SHM object
%@y:    current solution point
%@J:    Jacobian of the current solution point
%@DYN:  Dynamical System class object
%@AM:   Approximation Method subclass object
%
%@multipliers:      Lyapunov exponents
%@vectors:          empty so far (EQ: eigenvectors of Jacobian, PS: eigenvectors of monodromy matrix)
%@n_unstable:       number of unstable Lyapunov exponents
%@stability_flag:   flag indicating the successfull computation of the multipliers

function   [multipliers,vectors,n_unstable,stability_flag] = QPS_SHM_calc_stability(obj,y,J,DYN,AM)

    stability_flag = 1;
    dim = AM.n;
    nchar = AM.n_char;
    phi = AM.phi;
    n_auto = DYN.n_auto;
    
    try
    
        if(DYN.n_auto==0)
            Omega = DYN.non_auto_freq(y(end,1));
            if(Omega(1,1) > Omega(1,2));  index = [1,2]; else;  index = [2,1]; end              % Define integration variable to minimize integration time
        elseif(DYN.n_auto==1)
            index = [1,2];
            Omega(1,1) = DYN.non_auto_freq(y(end,1));
            Omega(1,2) = y(end-1,1);
        elseif(DYN.n_auto==2)
            index = [1,2];
            Omega(1,:) = y(end-2:end-1,1);
        end
        Ik = [0,2.*pi./Omega(1,index(1,1))];                                                     % Set time span for integration
        Theta = mod(phi(index(1,2),:) + Ik(2)*Omega(1,index(1,2)),2*pi);                         % Calculate end values of characteristics and map back to 0,2pi square
        [~,I(1,:)] = sort(Theta);                                                                % Sort remapped values in ascending order
    
        %%%%%%%%%%%%
        %%%%%%%%%%%%
        % dz_int_dz0 = obj.jacobi_int(y,Omega,AM,index);                                         % Ableitung dZ0~/dZ0
        %%%%%%%%%%%%
        %%%%%%%%%%%%
    
        if strcmpi(DYN.approx_method,'shooting')                                                % Solution has already been computed by shooting - no need to do it again
            for k=1:nchar                                                                       % Assemble Fundamental matricies
                r = I(1,k);
                %         FundMatrix(1:dim,1:dim,k) = J(dim*(k-1)+1:k*dim,dim*(r-1)+1:r*dim);%+dz_int_dz0(dim*(r-1)+1:r*dim,dim*(k-1)+1:k*dim);
                FundMatrix(1:dim,1:dim,r) = J(dim*(k-1)+1:k*dim,dim*(r-1)+1:r*dim);%+dz_int_dz0(dim*(r-1)+1:r*dim,dim*(k-1)+1:k*dim);
                FundMatrix(1:dim,1:dim,r) = J(dim*(k-1)+1:k*dim,dim*(r-1)+1:r*dim);
            end
    
            % Consider periodicity
            FundMatrixARRAY  = cat(3,FundMatrix,FundMatrix(:,:,1));                             % Add first fundamental matrix for periodicity
    
            for j=1:dim                                                                         % Function of the fundamental solution
                for i=1:dim
                    %             FundMatrixFUN{i,j} = griddedInterpolant([phi(index(1,2),:),2*pi],permute(FundMatrixARRAY(i,j,:),[1,3,2]),'spline');                      %Doppelte For-Schleife vektorisieren
                    FundMatrixFUN{i,j} = @(x) fnval(csape([phi(index(1,2),:),2*pi],permute(FundMatrixARRAY(i,j,:),[1,3,2]),'periodic'),x);
                end
            end
            % Modulate function as array output
            %     FundFUN = @(x)cellfun(@(c) c(x),FundMatrixFUN,'UniformOutput',false);
            FundFUN = @(x)cellfun(@(c) c(x),FundMatrixFUN,'UniformOutput',false);               % Make fundamental matrix function cell function
    
            % Calculate Ljapunov Exponents
            NOM = 2*10^4;                                                                       % Number of Multiplication
            x = 1:dim;
    
            FundNorm = eye(dim);                                                                % Take triangular matrix instead of identity matrix for columns in the jacobian with only one entry not equal to zero
    
            % All coordinates of mapping
            TCR = permute(mod(abs(((1:NOM)-1)*2*pi*(Omega(1,index(1,2)))/(Omega(1,index(1,1)))),2*pi),[1 3 2]);
            % All fundamental matrices of mapping
            FundARRAY = cell2mat(FundFUN(TCR));
            FundARRAYVal = zeros(size(FundARRAY));
    
            for i=1:NOM                                                                         % Mapping of perturbation
                FundARRAYVal(:,:,i) = FundARRAY(:,:,i)*FundNorm;                                % Multiplication: current fundamental soultion and normalized solution
                FundNorm = obj.GSOrthonormalization(FundARRAYVal(:,:,i));                       % Gram-Schmidt orthonormalisation
            end
    
            TimesFunc = @(x) FundARRAYVal(:,:,x)'*FundARRAYVal(:,:,x);                          % Multiplication for Gram determinant
            Vsqr = reshape(cell2mat(arrayfun(TimesFunc, 1:NOM,'UniformOutput',false)),dim,dim,[]);      % Ersetzen
    
            % Calculate volume
            detFunc = @(x,y) sqrt(det(Vsqr(1:x,1:x,y)));                                                    % Gram determinant
            Vol = reshape(arrayfun(detFunc,repmat(x,1,NOM),reshape(repmat(1:NOM,dim,1),1,[])),dim,[])';     % Volume of mapped perturbation
    
            % Catch errror
            if max(max(imag(Vol))) > sqrt(eps)
                error('!!! USERINFO: Time interval is too large for stability algorithm .... !!!')
            end
    
            LyExN = real(sum(log(Vol)))./(NOM*Ik(2));                                                       % Lyapunov exponents of n-th order
            Lyp_Values = sort([LyExN(1),diff(LyExN)],'descend');                                            % Calculate LEs of first order
    
        else %Recompute the solution with shooting algorithm to get the monodromy matrix
    
            error('Stability of Quasi-periodic solutions can only be calculated if Approximation method is "shooting"!')
    
        end
    
        % Analyse Lyapunov exponents
        [LEAbsSortVal,~] = sort(abs(Lyp_Values));                                                       % Sort LE by absolute value
        LEValwoZero = LEAbsSortVal(n_auto+1:end);                                                       % Extract zeros LE's due to autonomous frequencies
        LEVal = zeros(size(LEValwoZero,2),1);                                                           % Preallocation
    
    
        for i=1:size(LEValwoZero,2)                                                                     % Identify sign's of LE's
            Idx = [find(Lyp_Values == -LEValwoZero(i)),find(Lyp_Values == LEValwoZero(i))];
            LEVal(i) = Lyp_Values(Idx);                                                                 % Set LE's for stability determination
        end
    
        multipliers = LEVal;                                                                            % Set multipliers to calculated LE's
        [~,idx] = sort(abs(multipliers));                                                               % Ensure sorting according to norm
        multipliers = multipliers(idx);                                                                 % Sort LE's
    
        n_unstable = numel(find(multipliers>10*eps));                                                   % Find number of unstable multipliers
    
    catch   %Unlikely case for this computation
    
        stability_flag = 0;
    
    end
    
    vectors = [];
    
    [multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag);% Checks for NaN or Inf values
    

end
