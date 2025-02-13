% This function is a method of subclass AM_PS_FDM.
% It is needed for computing the stability via the shooting method and returns the state space vector z(theta=0).
% This vector is used as an initial condition by the shooting solver.
%
% @obj: ApproxMethod subclass AM_PS_FDM object
% @y:   Solution vector of the continuation in method space (method solution vector, autonomous frequencies, continuaton parameter)
% @DYN: DynamicalSystem class object
% @IC:  Initial condition vector in state space

function IC = getIC(obj,y,DYN,n_shoot)                                            

    dim = DYN.dim;                                          % Dimension of state space 
    n_int = obj.n_int;                                      % Number of intervals
    n_auto = DYN.n_auto;                                    % Number of autonomous frequencies

    if n_auto == 0
        omega = DYN.non_auto_freq(y(end,1));                % Get the frequency
    elseif n_auto == 1
        omega = y(end-1);                                   % Get the frequency
    end
    
    T = 2.*pi/omega;                                        % Periodic time
    t = linspace(0,T,n_int+1);                              % Time vector for interpolation

    Z = reshape(y(1:end-1-n_auto),dim,n_int);
    Z_interp = csape(t,[Z,Z(:,1)],'periodic');              % Spline interpolation of FD solution with periodic boundary conditions

    t_eval = linspace(0,T*(1-1/n_shoot),n_shoot);           % Time vector at shooting points
    Z_eval = fnval(Z_interp,t_eval);                        % Evaluate interpolation

    IC = reshape(Z_eval,dim*n_shoot,1);                     % The initial condition vector is the state space vector at the shooting points

end