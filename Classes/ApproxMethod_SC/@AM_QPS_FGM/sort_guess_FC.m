%Function is a method of ApproxMethod subclass AM_QPS_FGM and computes
%an initial value vector base on the opt_init option structure
% If there are less initial conditions than higher harmonics in obj.hmatrix,
% additional intial values are guessed
%
%@obj:  ApproxMethod class object of AM_PS_FGM
%@FC0:  Initial Fourier Coefficient vector
%
%@s:            update solution vector in approximation space (without auto-frequencies or continuation parameter)
%@hmatrix:      update matrix of higher harmonics

function [s,hmatrix] = sort_guess_FC(obj,DYN,FC0)
 

    tmp = reshape(FC0,DYN.dim,[]);
    c0 = tmp(:,1);
    cmatrix = tmp(:,2:0.5*(numel(tmp)/DYN.dim+1));          %Gatekeeper has ensured that this computation of dimension is always possible
    smatrix = tmp(:,0.5*(numel(tmp)/DYN.dim+1)+1:end);
 
    %Sort the hmatrix and the IVs accordingly (also account for case when there are more higher harmonics than Fourier coefficients)
    idx_tmp = min([size(cmatrix,2),size(obj.hmatrix,2)-1]);                 %-1 due to constant higher harmonic
    hhm = obj.hmatrix;
    
    hhm(:,find(prod(hhm==[0;0],1))) = [];                                                          %delete the zero frequency for now. Needed for sorting so that the indices of sorted hmatrix can be applied to cmatrix and smatrix
    [~,idx_tmp] = sort(vecnorm(hhm(:,1:idx_tmp),2,1));
    hh_tmp = hhm(:,idx_tmp);
    
    cmatrix = cmatrix(:,idx_tmp);
    smatrix = smatrix(:,idx_tmp);
    
    %Assure that the order iss [1;0] and then [0;1] (no real reason behind that): Since we sorted hmatrix already according to the Euclidean norm, these two frequency combinations will always be first in the matrix.
    [~,idx_tmp] = sort(hh_tmp(2,1:2));
    hh_tmp(:,1:2) = hh_tmp(:,idx_tmp);
    cmatrix(:,1:2) = cmatrix(:,idx_tmp);
    smatrix(:,1:2) = smatrix(:,idx_tmp);
    
    if size(cmatrix,2)==(obj.p_n_hh-1)
        s = [c0(:);cmatrix(:);smatrix(:)];         %For the case that obj.fc0 is given for all higher harmonics in hmatrix... this is a stupid way to do it. But for any other case it is easier to implement.
        hmatrix = [zeros(2,1),hh_tmp];
    else
        % We are guessing the additionally needed initial conditions based on the knowledge, that Fourier coefficients drop exponentially:
        % C_(n,m) = K/(n^b1*m^b2). K is set to be always C_(1,0) (that is deliberate).
    
        hh_guess = ((setdiff(hhm.',hh_tmp.','rows'))).';   %Get the higher harmonics which we want to guess.
        [~,idx_tmp] = sort(vecnorm(hh_guess,2,1));
        hh_guess = hh_guess(:,idx_tmp);
    
        %We've got four different zones for parameter generation and guessing:
        %  2D-Higher Harmonics
        %
        %          theta2
        %           |
        %   x x x x x x x x x
        %           |
        %   x x x x x x x x x
        %           |
        %           x x_x_x_x___ theta1
        %
        % a) Axis along theta1 hmatrix(n,0)
        % b) Axis along theta2 hmatrix(0,n)
        % c) Area in the first quadrant hmatrix(n>0,m>0)
        % d) Area in the second quadrant hmatrix(n<0,m>0)
        %
        % Ideally, we would determine and guess parameters for each of these
        % areas separately:
        % a) C_{n,0} = K_{1,0}/n^\beta_{a,1}
        % b) C_{0,n} = K_{0,1}/n^\beta_{b,1}
        % c) C_{n,m} = K_{c}/(n^\beta_{c,1}m^\beta_{c,2})
        % d) C_{n,m} = K_{d}/(n^\beta_{d,1}m^\beta_{d,2})
        %
        % But most likely, we will not have enough parameters to do so. So:
    
        %Let's separate into the 4 cases:
        hh_guess_a  = hh_guess(:,hh_guess(2,:)==0);
        hh_tmp_a    = hh_tmp(:,hh_tmp(2,:)==0);
        cmatrix_a   = cmatrix(:,hh_tmp(2,:)==0);
        smatrix_a   = smatrix(:,hh_tmp(2,:)==0);
    
        hh_guess_b  = hh_guess(:,hh_guess(1,:)==0);
        hh_tmp_b    = hh_tmp(:,hh_tmp(1,:)==0);
        cmatrix_b   = cmatrix(:,hh_tmp(1,:)==0);
        smatrix_b   = smatrix(:,hh_tmp(1,:)==0);
    
        hh_guess_c  = hh_guess(:,logical((hh_guess(1,:)>0).*(hh_guess(2,:)~=0)));
        hh_tmp_c    = hh_tmp(:,logical((hh_tmp(1,:)>0).*(hh_tmp(2,:)~=0)));
        cmatrix_c   = cmatrix(:,logical((hh_tmp(1,:)>0).*(hh_tmp(2,:)~=0)));
        smatrix_c   = smatrix(:,logical((hh_tmp(1,:)>0).*(hh_tmp(2,:)~=0)));
    
        hh_guess_d  = hh_guess(:,(hh_guess(1,:)<0));
        hh_tmp_d    = hh_tmp(:,(hh_tmp(1,:)<0));
        cmatrix_d   = cmatrix(:,(hh_tmp(1,:)<0));
        smatrix_d   = smatrix(:,(hh_tmp(1,:)<0));
    
        % These are the default values... may get overwritten
        C_a = abs(cmatrix_a(:,1));  %Always known, since Gatekeeper ensured that there is [1,0] and [0,1]
        S_a = abs(smatrix_a(:,1));  %Always known, since Gatekeeper ensured that there is [1,0] and [0,1]
        C_b = abs(cmatrix_b(:,1));  %Always known, since Gatekeeper ensured that there is [1,0] and [0,1]
        S_b = abs(smatrix_b(:,1));  %Always known, since Gatekeeper ensured that there is [1,0] and [0,1]
        C_c = abs(mean([cmatrix_a,cmatrix_b],2));
        S_c = abs(mean([smatrix_a,smatrix_b],2));
        C_d = C_c;
        S_d = S_c;
    
        [C_beta_a,S_beta_a,C_beta_b,S_beta_b,C_beta_cI,C_beta_cII,S_beta_cI,S_beta_cII,C_beta_dI,C_beta_dII,S_beta_dI,S_beta_dII] = deal(ones(size(C_a)));
    
        % Find the exponents for area A
        if size(hh_tmp_a,2)>1
            C_beta_a = guess_ab(cmatrix_a,hh_tmp_a);
            S_beta_a = guess_ab(smatrix_a,hh_tmp_a);
            if size(hh_tmp_b,2)==1  %Use the found values for area B if there is not enough data
                C_beta_b = C_beta_a; S_beta_b = S_beta_a;
            end
        end
        % Find the exponents for area B
        if size(hh_tmp_b,2)>1
            C_beta_b = guess_ab(cmatrix_b,hh_tmp_b);
            S_beta_b = guess_ab(smatrix_b,hh_tmp_b);
            if size(hh_tmp_a,2)==1 %Use the found values for area A if there is not enough data
                C_beta_a = C_beta_b; S_beta_a = S_beta_b;
            end
        end
        % Find the coefficient and exponents for area C
        if size(hh_tmp_c,2)>1
            [C_c,C_beta_cI,C_beta_cII] = guess_cd(cmatrix_c,hh_tmp_c);
            [S_c,S_beta_cI,S_beta_cII] = guess_cd(smatrix_c,hh_tmp_c);
            if size(hh_tmp_d,2)<3 %Use the found values for area D if there is not enough data
                C_beta_dII = C_beta_cII; S_beta_dII = S_beta_cII;
                if size(hh_tmp_d,2)<2
                    C_beta_dI = C_beta_cI; S_beta_dI = S_beta_cI;
                    if size(hh_tmp_d,2)<1
                        C_d = C_c; S_d = S_c;
                    end
                end
            end
        end
    
        % Find the coefficient and exponents for area D
        if size(hh_tmp_d,2)>1
            [C_d,C_beta_dI,C_beta_dII] = guess_cd(cmatrix_d,hh_tmp_d);
            [S_d,S_beta_dI,S_beta_dII] = guess_cd(smatrix_d,hh_tmp_d);
            if size(hh_tmp_c,2)<3 %Use the found values for area C if there is not enough data
                C_beta_cII = C_beta_dII; S_beta_cII = S_beta_dII;
                if size(hh_tmp_c,2)<2
                    C_beta_cI = C_beta_dI; S_beta_cI = S_beta_dI;
                    if size(hh_tmp_c,2)<1
                        C_c = C_d; S_c = S_d;
                    end
                end
            end
        end
    
        guess_cmatrix_a = C_a./hh_guess_a(1,:).^C_beta_a;
        guess_smatrix_a = S_a./hh_guess_a(1,:).^S_beta_a;
        guess_cmatrix_b = C_b./hh_guess_b(2,:).^C_beta_b;
        guess_smatrix_b = S_b./hh_guess_b(2,:).^S_beta_b;
    
        guess_cmatrix_c = C_c./(hh_guess_c(1,:).^C_beta_cI.*hh_guess_c(2,:).^C_beta_cII);
        guess_smatrix_c = S_c./(hh_guess_c(1,:).^S_beta_cI.*hh_guess_c(2,:).^S_beta_cII);
    
        guess_cmatrix_d = C_d./(abs(hh_guess_d(1,:)).^C_beta_dI.*abs(hh_guess_d(2,:)).^C_beta_dII);
        guess_smatrix_d = S_d./(abs(hh_guess_d(1,:)).^S_beta_dI.*abs(hh_guess_d(2,:)).^S_beta_dII);
    
    
    
        %Guess the additional initial conditions.
        tmp_cmatrix =  [cmatrix,guess_cmatrix_a,guess_cmatrix_b,guess_cmatrix_c,guess_cmatrix_d];
        tmp_smatrix =  [smatrix,guess_smatrix_a,guess_smatrix_b,guess_smatrix_c,guess_smatrix_d];
    
    
        %Sort the complete h-matrix again. Why? Maybe ferquencies needed to be guessed which are not higher than the hh_tmp frequencies but in between (unlikely case... but now accounted for).
        hhm2 = [hh_tmp,hh_guess_a,hh_guess_b,hh_guess_c,hh_guess_d];
        [~,idx_tmp] = sort(vecnorm(hhm2,2,1));
        hmatrix = [zeros(2,1),hhm2(:,idx_tmp)];
    
        cmatrix = tmp_cmatrix(:,idx_tmp);
        smatrix = tmp_smatrix(:,idx_tmp);
    
        s = [c0(:);cmatrix(:); smatrix(:)];
    
    end
    

    
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Function used for guessing the coefficient and exponents of the areas A and B
function beta = guess_ab(cmatrix,hh_tmp)
    
    %hh_tmp(:,end) should be the largest harmonic due to sorting. Max ensures, that the non-zero element is chosen. There are no negative elements in this matrix.
    beta = mylog(abs(cmatrix(:,1))./abs(cmatrix(:,end)),max(hh_tmp(:,end)));
    
    beta = max([ones(size(beta)),beta],[],2);   %This ensures that any erorrs or numerical failures will not result in growing Fourier coefficients.
    end
    
    %Function used for guessing the coefficient and exponents of the areas C and D
    function [C,beta1,beta2] = guess_cd(cmatrix,hh_tmp)
    
    if size(hh_tmp,2)== 1
        %beta1, beta2 are both by default 1, thus the constants are:
        C = abs(cmatrix(:,1));
        beta1 = ones(size(C));
        beta2 = beta;
    elseif size(hh_tmp,2)==2
        C0 = abs(cmatrix(:,1));
        beta10 = ones(size(C0));
    
        try
            options = optimoptions('fsolve','Display','off');
            x = fsolve(@(x)qp_beta1(x,cmatrix(:,1),hh_tmp(:,1),cmatrix(:,2),hh_tmp(:,2)),[beta10;C0],options);
        catch
            x = [beta10;C0];
            warning('S.th. went wront with predicting higher harmonics');
        end
        beta1 = x(1:0.5*(numel(x)));
        beta2 = ones(size(beta1));
        C = x(0.5*numel(x)+1:end);
    else
        if size(hh_tmp,2)>3
            idx = mydist(hh_tmp); %This give back the indices of the three elements farest apart from each other
        else
            idx = [1,2,3];
        end
        C0 = abs(cmatrix(:,1));
        beta10 = ones(size(C0));
        beta20 = ones(size(C0));
        try
            options = optimoptions('fsolve','Display','off');
            x = fsolve(@(x)qp_beta12(x,cmatrix(:,idx(1)),hh_tmp(:,idx(1)),cmatrix(:,idx(2)),hh_tmp(:,idx(2)),cmatrix(:,idx(3)),hh_tmp(:,idx(3))),[beta10;beta20;C0],options);
        catch
            x = [beta10;beta20;C0];
            warning('S.th. went wront with predicting higher harmonics');
        end
        beta1    =          x(1:1/3*numel(x),:);
        beta2    =          x(1/3*numel(x)+1:2/3*numel(x),:);
        C        =          x(2/3*numel(x)+1:end,:);
    
    end
    
    %Values smaller than 1 do not seem to be uncommon.
    beta1 = max([zeros(size(beta1)),beta1],[],2);   %This ensures that any erorrs or numerical failures will not result in growing Fourier coefficients.
    beta2 = max([zeros(size(beta1)),beta2],[],2);   %This ensures that any erorrs or numerical failures will not result in growing Fourier coefficients.
    beta1(beta1==0) = 1;
    beta2(beta2==0) = 1;
    
end
    
    
function x = mylog(y,p) %logarithm from y w.r.t. base p
    
    x = log2(y)./log2(p);
    
end
    
    %Function defining the non-linear problem for finding the constant and beta1
function r = qp_beta1(x,Knm,h1,Kpell,h2)
    
    beta1    =        x(1:0.5*numel(x),:);
    Ktilde  =         x(0.5*numel(x)+1:end,:);
    
    h1 = abs(h1);
    h2 = abs(h2);
    Knm = abs(Knm);
    Kpell = abs(Kpell);
    
    r1 = Knm.*h1(1).^beta1.*h1(2) - Ktilde;
    r2 = Kpell.*h2(1).^beta1.*h2(1) - Ktilde;
    
    r = [r1;r2];
    
end
    
    %Function defining the non-linear problem for finding the constant and beta1 and beta2
function r = qp_beta12(x,K1,h1,K2,h2,K3,h3)
    
    beta1    =        x(1:1/3*numel(x),:);
    beta2    =         x(1/3*numel(x)+1:2/3*numel(x),:);
    Ktilde   =         x(2/3*numel(x)+1:end,:);
    
    h1 = abs(h1);
    h2 = abs(h2);
    h3 = abs(h3);
    K1 = abs(K1);
    K2 = abs(K2);
    K3 = abs(K3);
    
    r1 = K1.*h1(1).^beta1.*h1(2).^beta2 - Ktilde;
    r2 = K2.*h2(1).^beta1.*h2(2).^beta2 - Ktilde;
    r3 = K3.*h3(1).^beta1.*h3(2).^beta2 - Ktilde;
    
    r = [r1;r2;r3];
end
    
    
function idx = mydist(hh_tmp)
    %Find the three elements, which are farest apart from each other.
    %Find two harmonics which have the maximum distance from each other.
    %I do not know if this is smart or not... maybe find a better solution.
    %TODO: Find smater way to do that
    
    dist = zeros(size(hh_tmp,2),size(hh_tmp,2),size(hh_tmp,2));
    for kk = 1:size(hh_tmp,2) %rows
        for jj = kk+1:size(hh_tmp,2)
            for ii = jj+1:size(hh_tmp,2)
                tmp(1,1) = norm(hh_tmp(:,kk)-hh_tmp(:,jj));
                tmp(2,1) = norm(hh_tmp(:,jj)-hh_tmp(:,ii));
                tmp(3,1) = norm(hh_tmp(:,kk)-hh_tmp(:,ii));
                dist(kk,jj,ii) = mean(tmp);
            end
        end
    end
    
    tmp = dist(:);
    [~,idx_tmp] = sort(tmp,'descend');
    c= find(dist==tmp(idx_tmp(1)));
    [idx(1) idx(2) idx(3)] = ind2sub(size(dist),c(1));  %c(1), if there are multiple distances of the same size.
    
    
end
    

