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

    z0 = DYN.opt_init.ic;               % Initial point in state space
    n_auto = DYN.n_auto;                % Number of autonomous frequencies
    n_shoot = obj.n_shoot;              % Number of shooting points
    

    if n_auto == 0                      % Non-autonomous system
    
        if numel(z0) < dim*n_shoot      % Numerical time integration necessary to obtain the missing shooting points
                    
            mu = param{DYN.act_param};                                                      % Continuation parameter
            T = 2*pi/DYN.non_auto_freq(mu);                                                 % Periodic time
            T_int = linspace(0,T*(1-1/n_shoot),n_shoot);                                    % Time vector of shooting points
                    
            [~,Z] = obj.solver_function(@(t,z) Fcn(t,z,param), T_int, z0, obj.odeOpts);     % Integration to get required shooting points

            if n_shoot == 2             % In this case, t is a [1x2] vector and is therefore interpreted as [t_start, t_end] by the solver
                Z = Z([1,end],:);       % Thus, Z stores the values at all computed time points, but we only need the Z values at t
            end
            obj.iv = reshape(Z.',dim*n_shoot,1);                                            % Reshape to a vector

        elseif numel(z0) == dim*n_shoot

            obj.iv = z0;                % z0 can be used directly as initial values for the shooting points

        end


    elseif n_auto == 1                  % Autonomous system

        omega = DYN.auto_freq;          % Autonomous frequency
        T = 2*pi./omega;                % Periodic time

        T_int = linspace(0,T*(1-1/(n_shoot)),n_shoot);                                      % Time vector for the integration

        if numel(z0) ~= dim*n_shoot     % Numerical time integration necessary to obtain the missing shooting points

            [~,Z] = obj.solver_function(@(t,z) Fcn(t,z,param), T_int, z0, obj.odeOpts);     % Integration to get required shooting points (and evaluation points for integral phase condition)

            if n_shoot == 2             % In this case, T_int is a [1x2] vector and is therefore interpreted as [t_start, t_end] by the solver
                Z = Z([1,end],:);       % Thus, Z stores the values at all computed time points, but we only need the Z values at T_int
            end

            obj.iv = [reshape(Z.',dim*n_shoot,1); omega];                                   % Build the initial value
    
        elseif numel(z0) == dim*n_shoot

            obj.iv = [z0; omega];       % z0 can be used directly as initial values for the shooting points
            
        end

    end

end