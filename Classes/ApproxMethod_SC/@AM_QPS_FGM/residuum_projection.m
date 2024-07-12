%This function gives back an estimation of the error in frequency space for the current 
%approximation to the solution
%
%@obj:  ApproxMethod type object
%@y:    Current solution vector of the continuation
%@DYN:  DynamicalSystem class object
%
%@proj: projected residuum 
function proj = residuum_projection(obj,y,DYN)

   %(for some reason preallocating these variable is way faster
    % than using them directly)

    n_fft = obj.n_fft;              %FFT resolution    
    p_n_hh = obj.p_n_hh;            %number of higher harmonics
    p_chf = obj.p_chf;              %complex argument for the fourier series evaluation
    hmatrix = obj.hmatrix;          %Matrix of higher harmonics

    n_auto = DYN.n_auto;            %Number of autonomous frequencies
    dim = double(DYN.dim);          %dimension of state space. conversion to double is needed for computation of the indices. 
    Fcn = DYN.rhs;                  %right-hand side of ODE

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

    %% Projection with FFT
    FC = [s(1:dim);s((dim+1):(p_n_hh)*dim)-1i.*s(((p_n_hh)*dim+1):end)];   %Assemble complex Fourrier coefficient vector
    FCtemp = reshape(FC,[dim,p_n_hh]);                                     %Reshape the vector for building  up the series
    FS = real(FCtemp*p_chf);                                               %Compute the real Fourier series with the argument from the class

    f = reshape((Fcn(p_arg_val(1:(2-n_auto),:),FS,param)).',[n_fft,n_fft,dim]);  %Evaluate the rhs and reshape it: Hypertime_1 x Hypertime_2 x nDoF. 
                                                                                    %p_arg_val(1:(2-n_auto)): In the mixed or pure autonomous case,
                                                                                    %no explicit hyper-time arguments are needed. Since the non_autonomous 
                                                                                    %frequencies are always first, this works
%     dZ = reshape((1i.*Gamma(1,1).*FCtemp*(hmatrix(1,:).'.*p_chf) + 1i.*Gamma(2,1).*FCtemp*(hmatrix(2,:).'.*p_chf)).',[n_fft,n_fft,dim]); 

    dZ = reshape((1i.*repmat(sum(repmat(Gamma,1,p_n_hh).*hmatrix,1),dim,1).*FCtemp*p_chf).',[n_fft,n_fft,dim]); %Compute the left hand side with the derivatives. Transposition is needed so that reshape works correctly

    proj = (fftshift(fftshift((fft2(real(dZ-f)))./(n_fft^2),1),2));               %Do the projection and shift

end