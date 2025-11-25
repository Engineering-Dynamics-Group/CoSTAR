% Method of SOL_QPS_SHM: This method (re-)calculates a quasi-periodic manifold of a shooting solution
%
% @obj:       Solution subclass object
% @DYN:       DynamicalSystem object
% @options:   options structure for postprocessing solutions
%
% @s_hypertime: Hypertime solution array: This must(!) be a [res_1 x res_2 x dim x n_evals] dimensional array !!!
% @mu:          Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @hypertime:   Array of the time points: This must(!) be a [res_1 x res_2 x 2 x n_evals] dimensional array !!!
% n_evals:      Number of curve points to be evaluated 

function [s_hypertime,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

%% Parameters
Fcn = DYN.rhs;                                                          % Set Fcn to RHS of ODE
dim = obj.n;                                                            % Get the state-space dimension
param = DYN.param;                                                      % Get the param-array
n_char = obj.n_char;                                                    % Get the number of characteristics
index = options.index;                                                  % Get the index of the solution
n_evals = numel(index);                                                 % Get the number of solution points to be calculated
mu = obj.mu(1,index);                                                   % Get the mu-values
s1 = reshape(obj.s(:,index),dim,n_char,n_evals);                        % Get solution vectors at desired mu-values and reshape them
freq = obj.freq(:,index);                                               % Get the corresponding frequencies
phi = obj.phi;                                                          % Get the phi-values
if size(options.resolution,2) == 2
    res = options.resolution;                                           % Desired resolution of output s if user defined [1x2] array
else
    res = repmat(options.resolution,1,2);                               % Desired resolution of output s if user defined scalar
end
res_t = max(res);                                                       % Set resolution for time integration


%% Initialization
s_hypertime = zeros(res(1),res(2),dim,n_evals);                         % Initialize interpolated values on grid

% Calculate the hypertime array for the output
theta_eval_1 = linspace(0, 2*pi, res(1));                               % Coordinates of the output in theta_1-direction
theta_eval_2 = linspace(0, 2*pi, res(2));                               % Coordinates of the output in theta_2-direction
[Theta_2, Theta_1] = meshgrid(theta_eval_2,theta_eval_1);               % Theta_2 is the first output of meshgrid so that hypertime(:,:,2) has the desired structure
hypertime(:,:,1) = Theta_1;                                             % Theta_1 values as a matrix. First row is the zero vector, second row only consists of values 2*pi/res etc.
hypertime(:,:,2) = Theta_2;                                             % Theta_2 values as a matrix. First column is the zero vector, second column only consists of values 2*pi/res etc.
hypertime = repmat(hypertime,[1,1,1,n_evals]);                          % The matrices are the same for each evaluation


%% Transform time solution to hyper-time
for k = 1:n_evals                                                       % For all curve-points given by index

    % Time integration
    param{DYN.act_param} = mu(k);                                       % Set active parameter
    Omega = freq(:,k).';                                                % Set frequencies
    
    PHI(1,:) = phi(1,:)./Omega(1,1);                                    % Calculate Spacing of characteristics
    PHI(2,:) = phi(2,:)./Omega(1,2);
    if(DYN.n_auto==0)                                                   % If n_auto=0 choose best integration "direction"
        if(Omega(1,1) > Omega(1,2));  index_Om = [1,2]; else;  index_Om = [2,1]; end
    elseif(DYN.n_auto==1)
        index_Om = [1,2];
    elseif(DYN.n_auto==2)
        index_Om = [1,2];  
    end
    Ik = [0,2.*pi./Omega(1,index_Om(1,1))];                             % Set integation interval
    
    tspan = linspace(Ik(1),Ik(2),res_t);                                % Set grid for time integration
    [time,s_time] = obj.solver_function(@(t,x)obj.FcnWrapper_SOL_ODE2(t,x,@(t,x)Fcn(t,x,param),PHI),tspan,s1(:,:,k),obj.odeOpts); % Do the time-intergation
    s_time = reshape(s_time,[res_t,dim,n_char]);                        % Reshape the solution values of the integration into a 3D array
    
    % Interpolation
    Iota_1 = Omega(1).*time;                                            % Get hyper-time coordinates theta1
    Iota_2 = Omega(2).*time;                                            % Get hyper-time coordinates theta2
    Iota_3 = mod(repmat(Iota_1,1,n_char)+phi(1,:),2*pi);                % Add the phase difference of the starting points of the charcteristics and remap values into square [0,2+pi]
    Iota_4 = mod(repmat(Iota_2,1,n_char)+phi(2,:),2*pi);                % Add the phase difference of the starting points of the charcteristics and remap values into square [0,2+pi]

    Iota_11 = reshape(Iota_3,[n_char*res_t,1]);                         % Reshape theta_1-values of solution values into column vector
    Iota_22 = reshape(Iota_4,[n_char*res_t,1]);                         % Reshape theta_2-values of solution values into column vector
    Iota_33 = reshape(permute(s_time,[1,3,2]),[n_char*res_t,dim]);      % Reshape the solution values into a matrix: The i-th column stores the values along all characteristics for the i-th system dimension
    
    for j = 1:dim
        I3 = Iota_33(:,j);                                              % Get data for interpolation
        F = scatteredInterpolant(Iota_11,Iota_22,I3,'natural');         % Interpolate scattered data from time integration
        s_hypertime(:,:,j,k) = F(Theta_1,Theta_2);                      % Evaluate the interpolated data
    end

    % Interpolate boundaries of other torus coordinate
    s_hypertime(1,:,:,k)   = 0.5.*(s_hypertime(2,:,:,k)+s_hypertime(end-1,:,:,k));      % Interpolate other boundary which is not smoth due to shooting never set a solution point on the boundary
    s_hypertime(end,:,:,k) = s_hypertime(1,:,:,k);                                      % Set values of opposite boundary equal, due to periodic boundary conditions
    s_hypertime(:,1,:,k)   = 0.5.*(s_hypertime(:,2,:,k)+s_hypertime(:,end-1,:,k));      % Interpolate other boundary which is not smoth due to shooting never set a solution point on the boundary
    s_hypertime(:,end,:,k) = s_hypertime(:,1,:,k);                                      % Set values of opposite boundary equal, due to periodic boundary conditions

end     

end