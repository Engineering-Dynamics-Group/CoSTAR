% This is a method of subclass AM_PS_FGM and is used for comuting the stability via the multiple shooting method
% and maybe other purposes
% Method for getting an initial point in state space on the periodic solution orbit
%
%@obj:  object of AM_PS_FGM subclass
%@y:    solution vector of the continuation in method space (method solution vector, autonomous frequencies, continuaton parameter)
%@DYN:  DynamicalSystem class object
%@IC:   initial conditions vector in state space (state space variable AND potentially autonomous frequency)

function IC = getIC(obj,y,DYN,n_shoot)

    dim = DYN.dim;                                          % Dimension of the state space
    n_auto = DYN.n_auto;                                    % Number of autonomous frequencies
    s = y(1:end-1-DYN.n_auto);                              % Get the Fourier-Coefficients
    n_hh = (size(s,1)/dim-1)/2+1;                           % Compute the number of higher harmonics

    if n_auto == 0
        omega = DYN.non_auto_freq(y(end,1));                % Get the frequency
    elseif n_auto == 1
        omega = y(end-1);                                   % Get the frequency
    end

    FC = [s(1:dim,1);s((dim+1):(n_hh)*dim,1)-1i.*s(((n_hh)*dim+1):end)];        % Assemble complex Fourier vector
    FC = reshape(FC,dim,n_hh);

    T = 2.*pi/omega;                                        % Periodic time
    t = linspace(0,T*(1-1/n_shoot),n_shoot);                % Time vector at shooting points
    chf =  exp(1i.*omega.*obj.hmatrix'.*t);                 % Evaluate complex Fourier series at shooting points

    s_out = real(pagemtimes(FC,chf));                       % Extract real Fourier coefficients
    IC = reshape(s_out,dim*n_shoot,1);                      % Extract state space values for multiple shooting

end