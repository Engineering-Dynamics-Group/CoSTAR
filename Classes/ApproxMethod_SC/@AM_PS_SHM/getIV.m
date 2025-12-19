% This function is a method of the subclass AM_PS_SHM.
% It creates x_IV, which is part of the initial value y_IV = [x_IV; mu_IV], for fsolve to find the initial solution.
% x_IV consists of the method solution vector s_IV and, in the autonomous case, of the autonomous frequency: x_IV = [s_IV; omega_IV]
%
% @obj: ApproxMethod subclass AM_PS_FDM object
% @DYN: DynamicalSystem class object

function obj = getIV(obj,DYN)
    
    % Parameter
    dim = DYN.dim;                      % Dimension of the system
    Fcn = DYN.rhs;                      % RHS of ODE
    param = DYN.param;                  % Parameter array

    s0 = DYN.opt_init.ic;               % Initial point in state space or method solution vector
    n_auto = DYN.n_auto;                % Number of autonomous frequencies
    n_shoot = obj.n_shoot;              % Number of shooting points
    

    if numel(s0) == dim*n_shoot         % The complete method solution vector is already supplied

        obj.iv = s0;                    % s0 can be used directly as initial values for the shooting points

        
    elseif numel(s0) == dim                                 % Only z(t=0) is supplied - time integration necessary to obtain the missing shooting points

        if n_auto == 0                                      % Non-autonomous system
            mu = param{DYN.act_param};                      % Continuation parameter
            T = 2*pi/DYN.non_auto_freq(mu);                 % Periodic time
        elseif n_auto == 1
            T = 2*pi/DYN.auto_freq;                         % Periodic time
        end
        T_int = linspace(0,T*(1-1/n_shoot),n_shoot);        % Time vector of shooting points

        [~,Z] = obj.solver_function(@(t,z) Fcn(t,z,param), T_int, s0, obj.odeOpts);     % Integration to get required shooting points

        if n_shoot == 2                                     % In this case, t is a [1x2] vector and is therefore interpreted as [t_start, t_end] by the solver
            Z = Z([1,end],:);                               % Thus, Z stores the values at all computed time points, but we only need the Z values at t
        end

        obj.iv = reshape(Z.',dim*n_shoot,1);                % Reshape to a vector

    
    else                                                                % All other cases: n_shoot_s0 = numel(s0)/dim - Interpolate s0

        n_shoot_s0 = numel(s0)/dim;                                     % Gatekeeper assures that n_shoot_s0 is an integer

        theta_interp = 0 : 2*pi/n_shoot_s0 : 2*pi;                      % theta values corresponding to the shooting points of s0
        Z_interp = [reshape(s0,dim,n_shoot_s0), s0(1:dim)];             % Reshape s0 to Z = [z(0), ... , z(2*pi-DeltaTheta)] and add z(0) (periodic condition)
        Z_csape = csape(theta_interp, Z_interp, 'periodic');            % Interpolation with periodic end conditions (1st and 2nd derivative). Output of csape is a struct (piecewise polynomial form)

        theta_eval = linspace(0,2*pi*(1-1/n_shoot),n_shoot);            % theta values for the evaluation (at the desired shooting points)
        Z_eval = fnval(Z_csape, theta_eval);                            % Evaluate the interpolated data
        obj.iv = reshape(Z_eval,dim*n_shoot,1);                         % Reshape to a vector

    end


    if n_auto == 1                              % If the system is autonomous

        obj.iv = [obj.iv; DYN.auto_freq];       % Add the autonomous frequency to the initial value

    end

end