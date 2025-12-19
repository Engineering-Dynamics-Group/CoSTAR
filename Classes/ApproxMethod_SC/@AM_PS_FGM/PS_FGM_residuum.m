%This function is a method of subclass AM_PS_FGM. This function computes
%the residuum with a Galerkin projection ones 
% for first order ODE systems
%
%@obj:  ApproxMethod subclass object
%@y:    solution curve vector (possibly containing auto frequency)
%@DYN:  DynamicalSystem class object
%
%@res:  residuum vector of evaluated ODE

function res = PS_FGM_residuum(obj,y,DYN)                            
        
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

    %% Allocate memory and declare persistent variables
    res = zeros(length(s)+DYN.n_auto,1);
    
    %% Projection with FFT

    FC = [s(1:dim);s((dim+1):(p_n_hh)*dim)-1i.*s(((p_n_hh)*dim+1):end)];    %Assemble complex Fourrier vector    
    FCtemp = reshape(FC,[dim,p_n_hh]);
    FS = real(FCtemp*p_chf);                                                %Real Fourier series

    proj = fft((Fcn(arg_val,FS,param)).').'./n_fft;

    Fnl = [real(proj(:,1)),2.*proj(:,(hmatrix(2:end)+1))];                 %Select the terms of the Fourier transform, which correspond to the hmatrix. "+1", since 0th frequency is the first entry in the vector.
    
    cmpl_res = 1i.*Gamma.*reshape(repmat(hmatrix,[dim,1]),[],1).*FC - Fnl(:);   %complex residuum equation

    %Assemble the residuum vector
    res(1:dim,1)                                            =   real(cmpl_res(1:dim));           %constant term
    res((dim+1):p_n_hh*dim,1)                               =   real(cmpl_res((dim+1):end,1));   %cosine term
    res((p_n_hh*dim+1):(2.*p_n_hh-1)*dim,1)                 =  -imag(cmpl_res((dim+1):end,1));   %sine term
   
    %% Add phase condition (only if autonomous frequencies are present)    
    if DYN.n_auto
        res(end,1) = obj.phase_condition(FCtemp,DYN);
    end

    %% Expand the residuum if system is conservative
    if isfield(DYN.system,'first_integral')
        I = DYN.system.first_integral;                              % Function of the first integral I = I(z)
        I_Z = I(FS,param);                                          % Evalute the first integral for all shooting points z_i
        IC = 1/n_fft*sum(I_Z) - param{end};                         % First Integral Constraint: I(s) = param{end} | Take the average of I_Z to get a single value for all shooting points
        res = [res; IC];                                            % Note: mean(I_Z) = 1/n_shoot*sum(I_Z), but the sum() function is somehow faster                                                      
    end


end