% This function is a method of subclass AM_QPS_FDM.
% It is needed for computing the stability via the shooting method and returns n_char_st state space vectors at z(0,theta_2).
% This vector is used as an initial condition by the shooting solver.
%
% @obj: ApproxMethod subclass AM_QPS_FDM object
% @y:   Solution vector of the continuation (method solution vector, autonomous frequencies, continuaton parameter)
% @DYN: DynamicalSystem class object
% @IC:  Initial condition vector in state space

function IC = getIC(obj,y,DYN,n_char_st)                                            

    % Parameters
    dim = DYN.dim;                                                      % Dimension of state space 
    n_auto = DYN.n_auto;                                                % Number of autonomous frequencies
    s = y(1:end-1-n_auto);                                              % Method solution vector
    mu = y(end);                                                        % Continuation parameter
    n_int_1 = obj.n_int_1;                                              % Number of intervals in theta_1-direction
    n_int_2 = obj.n_int_2;                                              % Number of intervals in theta_1-direction

    % Get the required state vectors of the solution
    Z = reshape(s,dim,n_int_1*n_int_2);                                 % Reshape the computed solution to a [dim x n_int_1*n_int_2] matrix
    if n_auto == 0                                                      % If system is fully non-autonomous
        Omega = DYN.non_auto_freq(mu);                                  % Get the frequencies
        if Omega(1) > Omega(2)                                          % If T1 < T2: Shooting integrates in theta_1 - direction
            Z_theta = Z(:,1:n_int_1:end);                               % Get the state vectors Z(theta_1=0,theta_2) as initial condition
            theta = linspace(0,2*pi,n_int_2+1);                         % theta_2 values of Z_theta (theta already includes 2*pi, whereas Z_theta does not yet!)
        else                                                            % If T2 <= T1: Shooting integrates in theta_2 - direction
            Z_theta = Z(:,1:n_int_1);                                   % Get the state vectors Z(theta_1,theta_2=0) as initial condition
            theta = linspace(0,2*pi,n_int_1+1);                         % theta_1 values of Z_theta (theta already includes 2*pi, whereas Z_theta does not yet!)
        end
    else                                                                % If system is (partly) autonomous: Shooting integrates in theta_1 - direction
        Z_theta = Z(:,1:n_int_1:end);                                   % Get the state vectors Z(theta_1=0,theta_2) as initial condition
        theta = linspace(0,2*pi,n_int_2+1);                             % theta_2 values of Z_theta (theta already includes 2*pi, whereas Z_theta does not yet!)
    end
    
    if size(Z_theta,2) == n_char_st                                     % No interpolation required - Z_theta values can be used directly
        Z_eval = Z_theta;
    else                                                                % Interpolation is required
        % Interpolation
        Z_interp = csape(theta,[Z_theta,Z_theta(:,1)],'periodic');      % Interpolate the state vectors (do not forget to append by the first element to have periodic boundary)
        % Evaluation
        theta_eval = linspace(0,2*pi*(1-1/(n_char_st+1)),n_char_st);    % theta values for the desired number of characteristics n_char_st
        Z_eval = fnval(Z_interp, theta_eval);                           % Evaluate interpolation at the desired theta values
    end
        
    % Initial states for n_char_st
    IC = reshape(Z_eval,dim*n_char_st,1);                               % Reshape the new initial states to a vector

end