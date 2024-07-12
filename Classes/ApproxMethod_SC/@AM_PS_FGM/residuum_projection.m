%This function gives back an estimation of the error in frequency space for the current 
%approximation to the solution
%
%@obj:  ApproxMethod object AM_PS_FGM
%@y:    Current solution vector of the continuation
%@DYN:  DynamicalSystem class object
%
% @proj: Projected Residuum

function proj = residuum_projection(obj,y,DYN)

   %(for some reason preallocating these variable is way faster
   % than using them directly)

    n_fft = obj.n_fft;      %FFT resolution
    dim = DYN.dim;          %dimension of state space
    p_n_hh = obj.p_n_hh;    %number of higher harmonics
    p_chf = obj.p_chf;      %complex argument for the fourier series evaluation
    hmatrix = obj.hmatrix;  %higher harmonics matrix
    Fcn = DYN.rhs;          %right-hand side of ODE

    s = y(1:(end-DYN.n_auto-1),1);           %cut off the autonomous frequencies and the curve parameter at the end of the constant vector (stacked Fourier coefficients)
    mu = y(end);

    %Scaling of the argument for the external excitation: cos(Omega*p_arg_val) --scale--> cos(Omega*arg_val)
    %This scaling is necessary, because the ODE function file is expected
    %to contain terms like cos(Omega*t) or something like that and we need
    %to account for that
    if DYN.n_auto
        arg_val = obj.p_arg_val./DYN.auto_freq;
        Gamma = y(end-DYN.n_auto);
    else
        arg_val = obj.p_arg_val./DYN.non_auto_freq(mu);
        Gamma = DYN.non_auto_freq(mu);
    end
    
    %Update the param array with the continuation parameter.
    param = DYN.param;
    param{DYN.act_param} = mu;

  
    %% Projection with FFT

    FC = [s(1:dim);s((dim+1):(p_n_hh)*dim)-1i.*s(((p_n_hh)*dim+1):end)];   %Assemble complex Fourrier vector    
    FCtemp = reshape(FC,[dim,p_n_hh]);
    
    FS = real(FCtemp*p_chf);                                              %Real Fourier series for z
    dZ = real(1i.*Gamma.*FCtemp*(hmatrix.'.*obj.p_chf));                  %Real Fourier series for dz/dtheta1

    proj = fft((dZ-Fcn(arg_val,FS,param)).').'./n_fft;                    %dz/dt - f(z,t) in frequency space


end