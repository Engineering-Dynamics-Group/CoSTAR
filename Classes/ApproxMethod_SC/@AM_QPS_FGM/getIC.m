% This is a method of subclass AM_QPS_FGM and is used for comuting the stability via the shooting method
% and maybe other purposes
% Method for getting an initial point in state space on the periodic solution orbit
%
% @obj:  object of AM_QPS_FGM subclass
% @y:    solution vector of the continuation (method solution vector, autonomous frequencies, continuaton parameter)
% @DYN:  DynamicalSystem class object
% @IC:   initial conditions vector in state space (state space variable AND potentially autonomous frequency)

function IC = getIC(obj,y,DYN,n_char_st)

    % Parameter
    dim = DYN.dim;                                                      % Dimension of the state space
    n_auto = DYN.n_auto;                                                % Number of autonomous frequencies
    s = y(1:end-1-n_auto);                                              % Get the Fourier-Coefficients
    mu = y(end);                                                        % Continuation parameter
    n_hh = (size(s,1)/dim-1)/2+1;                                       % Compute the number of (higher) harmonics
    hmatrix = obj.hmatrix;

    % Assemble complex Fourier vector
    FC = [s(1:dim); s((dim+1):(n_hh)*dim) - 1i.*s(((n_hh)*dim+1):end)];
    FC_mat = reshape(FC,dim,n_hh);

    % Get the required hypertime values
    theta = linspace(0,2*pi*(1-1/(n_char_st+1)),n_char_st);             % theta values for the desired number of characteristics n_char_st
    if n_auto == 0                                                      % If system is fully non-autonomous
        Omega = DYN.non_auto_freq(mu);                                  % Get the frequencies
        if Omega(1) > Omega(2)                                          % If T1 < T2: Shooting integrates in theta_1 - direction
            hypertime = [zeros(1,n_char_st);                            % theta_1 values for evaluation
                               theta       ];                           % theta_2 values for evaluation
        else                                                            % If T2 <= T1: Shooting integrates in theta_2 - direction
            hypertime = [      theta;                                   % theta_1 values for evaluation
                         zeros(1,n_char_st)];                           % theta_2 values for evaluation
        end
    else                                                                % If system is (partly) autonomous: Shooting integrates in theta_1 - direction
        hypertime = [zeros(1,n_char_st);                                % theta_1 values for evaluation
                           theta       ];                               % theta_2 values for evaluation
    end
    
    % Compute cos(H1*theta_1 + H2*theta_2) and sin(H1*theta_1 + H2*theta_2)
    chf = exp(1i.*hmatrix'*hypertime);      % complex harmonic functions: real(chf(2,:)) gives the cosine with the first higher harmonic defined ...
                                            % in hmatrix (first one is 0). imag(chf(2,:)) correspondingly gives the same for the sine function

    % Compute the Fourier series
    IC = reshape(real(FC_mat*chf),dim*n_char_st,1);                     % size(FC_mat*chf) =  [dim x n_hh] * [n_hh x n_char_st] = [dim x n_char_st]

end