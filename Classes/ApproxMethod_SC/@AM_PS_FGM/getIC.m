% This is a method of subclass AM_PS_FGM and is used for comuting the stability via the shooting method
% and maybe other purposes
% Method for getting an initial point in state space on the periodic solution orbit
%
%@obj:  object of AM_PS_FGM subclass
%@y:    solution vector of the continuation in method space (method solution vector, autonomous frequencies, continuaton parameter)
%@DYN:  DynamicalSystem class object
%@IC:   initial conditions vector in state space (state space variable AND potentially autonomous frequency)

function IC = getIC(obj,y,DYN)                                            

    dim = DYN.dim;                                  %dimension of the state space 
    s = y(1:end-1-DYN.n_auto);                      %Get the Fourier-Coefficients
    n_hh = (size(s,1)/dim-1)/2+1;                   %Compute the number of higher harmonics
    
    %For t = 0, only the cosine term coefficients are relevant: z = C_0 + \sum_{k = 1}^{n_hh} (C_k \cos(k \Omega t) +  S_k \sin(k \Omega t)) = C_0 + \sum_{k = 1}^{n_hh} C_k  
    IC = sum(reshape([s(1:dim,1);s((dim+1):(n_hh)*dim,1)],dim,n_hh),2);      %Assemble complex Fourrier vector
 
end