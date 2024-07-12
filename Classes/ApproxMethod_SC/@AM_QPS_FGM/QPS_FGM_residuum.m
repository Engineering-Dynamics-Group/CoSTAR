%This function is a method of subclass AM_QPS_FGM. This function computes
%the residuum with a Galerkin projection ones 
% for first order ODE systems
%
%@obj:  ApproxMethod subclass object
%@y:    solution curve vector (possibly containing auto frequency)
%@DYN:  DynamicalSystem class object
%
%@res:  residuum vector of evaluated ODE

function res = QPS_FGM_residuum(obj,y,DYN)                            
        
    %(for some reason preallocating these variable is way faster
    % than using them directly)

    n_fft = obj.n_fft;              %FFT resolution    
    p_n_hh = obj.p_n_hh;            %number of higher harmonics
    p_chf = obj.p_chf;              %complex argument for the fourier series evaluation
    hmatrix = obj.hmatrix;          %higher harmonics matrix

    n_auto = DYN.n_auto;            %Number of autonomous frequencies
    dim = double(DYN.dim);          %dimension of state space. conversion to double is needed for computation of the indices. 
    Fcn = DYN.rhs;                  %right-hand side of ODE

    n0 = n_fft/2+1;                 %Number important for finding the correct higher harmonics in the projected function 
    s = y(1:(end-DYN.n_auto-1),1);  %cut off the autonomous frequencies and the curve parameter at the end of the constant vector
    mu = y(end);                    %Current continuation parameter
  

    %This scaling is necessary, because the ODE function file is expected
    %to contain terms like cos(Omega*t) or something like that and we need
    %to account for that

    if DYN.n_auto == 0  %both frequencies are non-autonomous
        Gamma = reshape(DYN.non_auto_freq(mu),2,1);     %Make sure it is a column vector
    elseif DYN.n_auto == 1 %one is autonomous, one is non-autonomous
        Gamma = [DYN.non_auto_freq(mu);y(end-DYN.n_auto:end-1)];
    else %DYN.n_auto == 2 both frequencies are autonomous
        Gamma = y(end-DYN.n_auto:end-1);
    end
    p_arg_val = obj.p_arg_val./Gamma;

    %Update the param array with the continuation parameter.
    param = DYN.param;
    param{DYN.act_param} = mu;

    %% Allocate memory and declare persistent variables
    res = zeros(length(s)+DYN.n_auto,1);
    
    %% Projection with FFT
    FC = [s(1:dim);s((dim+1):(p_n_hh)*dim)-1i.*s(((p_n_hh)*dim+1):end)];   %Assemble complex Fourrier coefficient vector
    FCtemp = reshape(FC,[dim,p_n_hh]);                                     %Reshape the vector for building  up the series
    FS = real(FCtemp*p_chf);                                               %Compute the real Fourier series with the argument from the class

    temp = reshape((Fcn(p_arg_val(1:(2-n_auto),:),FS,param)).',[n_fft,n_fft,dim]);  %Evaluate the rhs and reshape it: Hypertime_1 x Hypertime_2 x nDoF. 
                                                                                    %p_arg_val(1:(2-n_auto)): In the mixed or pure autonomous case,
                                                                                    %no explicit hyper-time arguments are needed. Since the non_autonomous 
                                                                                    %frequencies are always first, this works
    
    proj = (fftshift(fftshift((fft2(temp))./(n_fft^2),1),2));               %Do the projection and shift


    indices =  ((n0-1)*n_fft+n0).*ones(dim,p_n_hh) +...
        repmat(sum([hmatrix(1,:);hmatrix(2,:).*n_fft],1),dim,1) +...
        repmat((0:(dim-1))'.*n_fft^2,1,p_n_hh);                            %This expression generates the (linear) indices based on Ktmp to get the correct components
                                                                           %The indices are ordered in a matrix to provoke matrix return of components
                                                                           %First summand adds a constant number to all other, since the fftshift function is used and the constant
                                                                           %coefficient is now in the middle of the spectrum.
                                                                           %Second summand generates the indices of the first layer of the 3D matrix.     

    Fnl = proj(indices).*repmat(logical(sum(abs(hmatrix),1))+1,dim,1);     %proj(indices) gives the components; repmat(...) 
                                                                           % defines a matrix multiplying ervery component with 2
                                                                           %except the constant one.
    
    hmat_gamma = reshape(repmat(sum(repmat(Gamma,1,p_n_hh).*hmatrix,1),dim,1),[],1);    %scalar product between hmatrix row and gamma ready for multiplication with the fourier coeff vector
    cmpl_res = 1i.*hmat_gamma.*FC - Fnl(:);                                 %Build the complex residuum

    
    res(1:dim,1)                                =   real(cmpl_res(1:dim));          %constant term
    res((dim+1):p_n_hh*dim,1)                   =   real(cmpl_res((dim+1):end,1));  %cosine term
    res((p_n_hh*dim+1):(2.*p_n_hh-1)*dim,1)     =  -imag(cmpl_res((dim+1):end,1));  %sine term
    
 
%     %% Add phase condition (only if autonomous frequencies are present)    
    if n_auto
        res((end-DYN.n_auto+1):end,1) = obj.phase_condition(FCtemp,DYN);
    end


end