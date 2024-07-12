% Method of SOL_QPS_SHM: This method (re-)calculates the time-solution of solutions calculated with quasi-periodic-shooting algorithm
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s:       Time solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @t:       Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function  [s,mu,t] = evalsol_time(obj,DYN,options)

index = options.index;
N = numel(index);

for k =1:N

    idx     = index(k);                                                                                                                     % Get index
    freq    = obj.freq(:,idx);                                                                                                              % Get frequencies
    mu_temp =  obj.mu(1,idx);                                                                                                               % Get bifurcation parameter
    s_temp  = reshape(obj.s,[DYN.dim,obj.n_char,size(obj.mu,2)]);                                                                           % Reshape vector s 
    
    Fcn = DYN.rhs;                                                                                                                          % RHS of ODE
    param = DYN.param;                                                                                                                      % Set parameters
    param{DYN.act_param} = mu_temp;                                                                                                         % Set active parameter

    if isfield(options,'interval')          %If an integration interval was supplied by the user... use this interval - else: used the standard interval
        s0init = s_temp(:,1,idx);                                                                                                           % Set initial solution for time integration
        if(options.interval(1)~=0)                                                                                                          % If starting point of interval is not 0
            [~,s_init] = obj.solver_function(@(t,z)Fcn(t+[0;0],z,param),[0,options.interval(1)],s0init,obj.odeOpts);                        % Integrate up to starting point of interval
            s0 = s_init(end,:).';                                                                                                           % Set initial condition according to starting point
        else
            s0 = s_temp(:,1,idx);                                                                                                           % If starting point of interval is zero skip initial integration
        end
        tspan = linspace(options.interval(1),options.interval(2),options.resolution);                                                       % Set integration interval
    else
        tspan = linspace(0,2*pi,options.resolution);                                                             % Set integration interval
        s0 = s_temp(:,1,idx);
    end
    [t(:,1,k) ,s(:,:,k) ] = obj.solver_function(@(t,z)Fcn(t+[0;0],z,param),tspan,s0,obj.odeOpts);                                           % Time integration for given interval
end

%Get the mu values
mu = obj.mu(options.index);    %options.index is unique due to S.solget_up_index

end






