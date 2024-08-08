% This is a method of subclass AM_PS_FGM and is used for comuting the stability via the multiple shooting method
% and maybe other purposes
% Method for getting an initial point in state space on the periodic solution orbit
%
%@obj:  object of AM_PS_FGM subclass
%@y:    solution vector of the continuation in method space (method solution vector, autonomous frequencies, continuaton parameter)
%@DYN:  DynamicalSystem class object
%@IC:   initial conditions vector in state space (state space variable AND potentially autonomous frequency)

function IC = getIC(obj,y,DYN,n_shoot)

dim = DYN.dim;                                                              % Dimension of the state space
s = y(1:end-1-DYN.n_auto);                                                  % Get the Fourier-Coefficients
n_hh = (size(s,1)/dim-1)/2+1;                                               % Compute the number of higher harmonics
Omega = DYN.non_auto_freq(y(end,1));

FC = [s(1:dim,1);s((dim+1):(n_hh)*dim,1)-1i.*s(((n_hh)*dim+1):end)];        % Assemble complex Fourrier vector
FC = reshape(FC,dim,n_hh);

T = linspace(0,2*pi/Omega,n_shoot+1);                                       % Define shooting points
chf =  exp(1i.*Omega.*obj.hmatrix'.*T);                                     % Evaluate compex Fourier series at shooting points

s_out = real(pagemtimes(FC,chf));                                           % Extract real Fourier coefficients
IC = reshape(s_out(:,1:end-1),[dim*n_shoot, 1]);                            % Extract intitial values for multiple shooting

end