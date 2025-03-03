% This is a method of the Stability subclass ST_QPS_SHM
% It calculates the Lyapunov exponents and indicates if the solution is stable or unstable
%
% @obj:  ST_QPS_SHM object
% @y:    current solution point
% @J:    Jacobian of the current solution point
% @DYN:  Dynamical System class object
% @AM:   Approximation Method subclass object
%
% @multipliers:      Lyapunov exponents
% @vectors:          empty so far (EQ: eigenvectors of Jacobian, PS: eigenvectors of monodromy matrix)
% @n_unstable:       number of unstable Lyapunov exponents
% @stability_flag:   flag indicating the successfull computation of the multipliers

function   [multipliers,vectors,n_unstable,stability_flag] = QPS_SHM_calc_stability(obj,y,J,DYN,AM)


    %% Parameters
    stability_flag = 1;                             % Initialise stability flag
    n_auto = DYN.n_auto;                            % Number of autonomous frequencies
    dim = DYN.dim;                                  % Dimension of state space
    n_map = obj.n_map;                              % Number of mappings to calculate the Ljapunov exponents
    n_char_st = obj.n_char_st;                      % Number of characteristics used to obtain the mapping matrices
    AM_ST = obj.AM_ST;                              % Object of AM_QPS_SHM used to reshoot the solution



    %% We need to recompute the solution with shooting algorithm if shooting was not used OR shooting was used but n_char ~= n_char_st
    % This is done to get the mapping matrices for the desired number of characteristics n_char_st
    if ~(strcmpi(DYN.approx_method,'shooting') && (AM.n_char == n_char_st))

        % We are not solving the problem with a subspace constraint due to the following reasons:
        % 1) Any other constraint except for the natural condition (y(end)-mu0) would no guarantee meeting the value of mu0.
        % 2) Even the natural constraint does not exactly guarantee that, since it is only solved approximately. There are cases
        %    where fsolve aborts the solution process even without exactly meeting the natural constraint. This occurs when
        %    the difference y(end)-mu0 becomes very small, which is likely when iterating a fold bifurcation point.

        % Get an initial value for fsolve
        x0 = AM.getIC(y,DYN,n_char_st);             % Get the initial states from a solution computed by finite differences or Fourier-Galerkin
        if n_auto > 0
            x0 = [x0; y(end-n_auto:end-1)];         % Add the autonomous frequency(s) to the initial states
        end

        % Reshoot the solution
        mu = y(end);                                % Continuation parameter
        phi = AM_ST.phi;
        AM_ST.IF_up_res_data([x0;mu],DYN);          % Update properties for phase condition
        J = speye(dim*n_char_st+n_auto);            % Preallocate Jacobian
        try
            [x,~,stability_flag,~,J] = fsolve(@(x) obj.QPS_SHM_ST_residuum(x,mu), x0, obj.fsolve_opts);      % Reshoot the solution to get J
            y = [x; mu];
            AM = AM_ST;
        catch
            stability_flag = 0;                     % If shooting did fail for some reasons, flag is set to 0
        end

    
    % If solution was already computed by shooting method with n_char = n_char_st
    else

        phi = AM.phi;

    end



    %% Calculate Ljapunov exponents
    
    if stability_flag > 0

        try

            if DYN.n_auto==0
                Omega = DYN.non_auto_freq(y(end,1));
                if(Omega(1,1) > Omega(1,2));  index = [1,2]; else;  index = [2,1]; end          % Define integration variable to minimize integration time
            elseif DYN.n_auto==1
                index = [1,2];
                Omega(1,1) = DYN.non_auto_freq(y(end,1));
                Omega(1,2) = y(end-1,1);
            elseif DYN.n_auto==2
                index = [1,2];
                Omega(1,:) = y(end-2:end-1,1);
            end
            Ik = [0,2.*pi./Omega(1,index(1,1))];                                                % Set time span for integration
            Theta = mod(phi(index(1,2),:) + Ik(2)*Omega(1,index(1,2)),2*pi);                    % Calculate end values of characteristics and map back to 0,2pi square
            [~,I(1,:)] = sort(Theta);                                                           % Sort remapped values in ascending order

            %%%%%%%%%%%%
            dz_int_dz0 = obj.jacobi_int(y,Omega,AM,index);                                      % Derivative dZ0~/dZ0
            % disp(max(abs(dz_int_dz0),[],'all'))
            %%%%%%%%%%%%

            FundMatrix = zeros(dim,dim,n_char_st);                                              % Preallocation for loop
            for k=1:n_char_st                                                                   % Assemble Fundamental matricies
                r = I(1,k);
                FundMatrix(1:dim,1:dim,r) = J(dim*(k-1)+1:k*dim,dim*(r-1)+1:r*dim) + dz_int_dz0(dim*(r-1)+1:r*dim,dim*(k-1)+1:k*dim);
            end
            FundMatrixARRAY  = cat(3,FundMatrix,FundMatrix(:,:,1));                             % Add first fundamental matrix for periodicity

            FundMatrixFUN = cell(dim,dim);                                                      % Preallocation for loop
            for j=1:dim                                                                         % Function of the fundamental solution
                for i=1:dim
                    % FundMatrixFUN{i,j} = griddedInterpolant([phi(index(1,2),:),2*pi],permute(FundMatrixARRAY(i,j,:),[1,3,2]),'spline');                      %Doppelte For-Schleife vektorisieren
                    FundMatrixFUN{i,j} = @(x) fnval(csape([phi(index(1,2),:),2*pi],permute(FundMatrixARRAY(i,j,:),[1,3,2]),'periodic'),x);
                end
            end
            % Modulate function as array output
            % FundFUN = @(x)cellfun(@(c) c(x),FundMatrixFUN,'UniformOutput',false);
            FundFUN = @(x)cellfun(@(c) c(x),FundMatrixFUN,'UniformOutput',false);               % Make fundamental matrix function cell function

            % Calculate Ljapunov Exponents
            d = 1:dim;
            FundNorm = eye(dim);                                                                % Take triangular matrix instead of identity matrix for columns in the jacobian with only one entry not equal to zero

            % All coordinates of mapping
            TCR = permute(mod(abs(((1:n_map)-1)*2*pi*(Omega(1,index(1,2)))/(Omega(1,index(1,1)))),2*pi),[1 3 2]);
            % All fundamental matrices of mapping
            FundARRAY = cell2mat(FundFUN(TCR));
            FundARRAYVal = zeros(size(FundARRAY));

            for i=1:n_map                                                                       % Mapping of perturbation
                FundARRAYVal(:,:,i) = FundARRAY(:,:,i)*FundNorm;                                % Multiplication: current fundamental soultion and normalized solution
                FundNorm = obj.GSOrthonormalization(FundARRAYVal(:,:,i));                       % Gram-Schmidt orthonormalisation
            end

            TimesFunc = @(x) FundARRAYVal(:,:,x)'*FundARRAYVal(:,:,x);                          % Multiplication for Gram determinant
            Vsqr = reshape(cell2mat(arrayfun(TimesFunc, 1:n_map,'UniformOutput',false)),dim,dim,[]);            % Replace

            % Calculate volume
            detFunc = @(x,y) sqrt(det(Vsqr(1:x,1:x,y)));                                        % Gram determinant
            Vol = reshape(arrayfun(detFunc,repmat(d,1,n_map),reshape(repmat(1:n_map,dim,1),1,[])),dim,[])';     % Volume of mapped perturbation

            % Catch errror
            if max(max(imag(Vol))) > sqrt(eps)
                error('Time interval is too large for stability computation algorithm!')
            end

            LyExN = real(sum(log(Vol)))./(n_map*Ik(2));                                         % Lyapunov exponents of n-th order
            Lyp_Values = sort([LyExN(1),diff(LyExN)],'descend');                                % Calculate LEs of first order


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


        catch   % Unlikely case for this computation

            stability_flag = 0;
            multipliers = NaN(dim-n_auto,1);
            n_unstable = NaN;

        end

    else

        multipliers = NaN(dim-n_auto,1);
        n_unstable = NaN;

    end
    

    vectors = [];
    
    [multipliers,vectors,n_unstable] = obj.check_stability_values(multipliers,vectors,n_unstable,stability_flag);   % Checks for NaN or Inf values
    

end