% Method of SOL_PS_MSHM: This method (re-)calculates a periodic orbit of a multiple shooting solution
% starting from the iterated points on the periodic orbit
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s:       Time solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @t:       Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function  [S,mu,t] = evalsol_time(obj,DYN,options)

    % Parameter
    dim = DYN.dim;              % Dimension of the system
    Fcn = DYN.rhs;              % RHS of ODE
    param = DYN.param;          % Parameter array

    index = options.index;      % Solution index (options.index is unique due to S.solget_up_index)
    n_eval = numel(index);      % Number of solutions to evaluate
    res = options.resolution;   % Resolution


    % Initialization
    S = zeros(res,dim,n_eval);
    t = zeros(res,1,n_eval);
    mu = obj.mu(index);         % Get the mu values


    % Evaluate for all requested solutions
    for i = 1:numel(index)
        
        % Parameter
        freq = obj.freq(1,index(i));      % Angular frequency
        s = obj.s(:,index(i));            % Method solution vector
        param{DYN.act_param} = mu(i);     % Update parameter array
        n_shoot = numel(s)/dim;           % Number of shooting points

        T = 2*pi/freq;                                                          % Periodic time
        if isfield(options,'interval')                                          % If an integration interval was supplied by the user
            t(:,1,i) = linspace(options.interval(1),options.interval(2),res);   % Use the supplied interval
        else                                                                    % If an integration interval was not supplied by the user
            t(:,1,i) = linspace(0,T,res);                                       % Integration interval is one period
        end


        % Get the data for one period by reshooting the solution using the integration from multiple shooting
        dT = T/n_shoot;                                     % Time span for each shooting operation (integration)
        T_int = [0:dT:(n_shoot-1)*dT; dT:dT:n_shoot*dT].';  % Define time intervals for the shooting operation (integration)
        z0_mat = reshape(s,dim,n_shoot);                    % Reshape s to a matrix of size [dim x n_shoot]
        T_ODE = cell(n_shoot,1);                            % Initialize a cell array to store all time values of the integration
        Z_ODE = cell(n_shoot,1);                            % Initialize a cell array to store all solution values of the integration
        numel_T = zeros(n_shoot,1);                         % Stores the number of time values of the integration
        for k = 1:n_shoot
            [t_ode,z_ode] = obj.solver_function(@(t,z) Fcn(t,z,param), T_int(k,:), z0_mat(:,k), obj.odeOpts);
            T_ODE{k} = t_ode;                               % Save the time values
            Z_ODE{k} = z_ode;                               % Save the solution values
            numel_T(k) = numel(t_ode);                      % Save the number of time values
        end


        % Interpolate the data for one period
        % We first need to assemble a time (row) vector t_interp and a [dim x numel(t_interp)] matrix Z_interp of the solution values for the interpolation
        % This could be done already in the loop above by appending an existing vector/matrix with the new ODE data
        % However, we do not know the sizes of t_ode and z_ode beforehand, thus we cannot initialize the needed vector/matrix to assemble
        % This is unfavourable, because extending the size of an existing variable in a loop is bad coding and can be slow
        t_interp = zeros(1,sum(numel_T)-n_shoot+1);         % Initialize vector to store the ODE time values for one period
        Z_interp = zeros(dim,sum(numel_T)-n_shoot+1);       % Initialize matrix to store the ODE solution values for one period
        count = 0;                                          % Counter telling how many elements have already been assembled to t_interp
        for k = 1:n_shoot                                   % This loop is to assemble t_interp and Z_interp 
            t_interp((count+1):(count+numel_T(k)-1)) = T_ODE{k}(1:end-1);       % We do not need the last entry since t_ode{k}(end) = t_ode{k+1}(1)
            Z_interp(:,(count+1):(count+numel_T(k)-1)) = Z_ODE{k}(1:end-1,:)';  % Again we do not need the last entry
            count = count + numel_T(k) - 1;                 % Add the number of added elements to the counter
        end
        t_interp(end) = T_ODE{end}(end);                    % Get the last entry at t = T (not included by the for loop above)
        Z_interp(:,end) = Z_ODE{end}(end,:)';               % Get the last entry at t = T (not included by the for loop above)
        Z_csape = csape(t_interp,Z_interp,'periodic');      % Spline interpolation of the solution with periodic boundary conditions


        % Evaluate the interpolated data
        t_eval = mod(t(:,1,i), T);                          % Reduce the requested time values to the interval [0,T) (since the interpolation was done in [0,T])
        S(:,:,i) = fnval(Z_csape,t_eval)';                  % Evaluate interpolation

        
    end


end