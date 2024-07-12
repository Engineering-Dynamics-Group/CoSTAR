% Method of SOL_QPS_SHM: This method (re-)calculates a quasi-periodic manifold of a shooting solution
%
% @obj:       Solution subclass object
% @DYN:       DynamicalSystem object
% @options:   options structure for postprocessing solutions
%
% @s_hypertime: Hypertime solution array: This must(!) be a [options.resolution x options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:          Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @hypertime:   Array of the time points: This must(!) be a [options.resolution x options.resolution x 2 x n_evals] dimensional array !!!
% n_evals:      Number of curve points to be evaluated 

function [s_hypertime,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

%% Initialization
index = options.index;                                                                                                  % Get the index of the solution
reso = options.resolution;                                                                                              % Get the resolution
n_char = obj.n_char;                                                                                                    % Get the number of characteristics
dim = obj.n;                                                                                                            % Get the state-space dimension
n_idx = numel(index);                                                                                                   % Get the number of solution points to be calculated
num_mu = numel(obj.mu);                                                                                                 % Get the number of curve points
tmp = linspace(0,2*pi,reso);                                                                                            % Set the interval for the Hyper-Time
[xq,yq] = meshgrid(linspace(0,2*pi,reso));                                                                              % Make a grid according to resolution
n_mu = size(options.index,2);                                                                                           % Get number solutions to calculate
[Iota_3,Iota_4] = deal(zeros(reso,obj.n_char,n_mu));                                                                    % Initialize matrices for remapped values
[Iota_11,Iota_22] = deal(zeros(obj.n_char *options.resolution,n_mu));                                                   % Initialize reshaped vales
Iota_33 = zeros(obj.n_char *reso,obj.n,n_mu);                                                                           % Initialize reshaped vales
S = zeros(reso,reso,obj.n,n_mu);                                                                                        % Initialize interpolated values on grid
time = zeros(reso,n_char,num_mu);                                                                                       % Initialize reshaped time vector
s0 = zeros(reso,dim,n_char,num_mu);                                                                                     % Initialize reshaped solution vector

s1 = reshape(obj.s,dim,n_char,num_mu);                                                                                  % Reshape solution vector
Fcn = DYN.rhs;                                                                                                          % Set Fcn to RHS of ODE

%% Transform time solution to hyper-time
counter = 0;                                                                                                            % Set counter
for k = index                                                                                                           % For all curve-points given by index
    counter = counter + 1;                                                                                              % Add counter
    %%%%%%%%% Time integration
    param = DYN.param;                                                                                                  % Fetch parameters
    param{DYN.act_param} = obj.mu(1,k);                                                                                 % Set active parameter
    Omega = obj.freq(:,k).';                                                                                            % Set frequencies
    
    PHI(1,:) = obj.phi(1,:)./Omega(1,1);                                                                                % Calculate Spacing of characteristics
    PHI(2,:) = obj.phi(2,:)./Omega(1,2);
    if(DYN.n_auto==0)                                                                                                   % If n_auto=0 choose best integration "direction"
        if(Omega(1,1) > Omega(1,2));  index_Om = [1,2]; else;  index_Om = [2,1]; end
    elseif(DYN.n_auto==1)
        index_Om = [1,2];
    elseif(DYN.n_auto==2)
        index_Om = [1,2];  
    end
    Ik = [0,2.*pi./Omega(1,index_Om(1,1))];                                                                             % Set integation interval
    
    tspan = linspace(Ik(1),Ik(2),reso);                                                                                 % Set grid for time integration
    [t_temp,s_temp] = obj.solver_function(@(t,x)obj.FcnWrapper_SOL_ODE2(t,x,@(t,x)Fcn(t,x,param),PHI),tspan,s1(:,:,k),obj.odeOpts); % Do the time-intergation
    time(:,:,k) = repmat(t_temp,1,n_char);                                                                              % Reshape the time vector
    s0(:,:,:,k) = reshape(s_temp,[reso,dim,n_char]);                                                                    % Reshape the solution vector
    
    Iota_1 = obj.freq(1,k).*time(:,1,k);                                                                                % Get hyper-time coordinates theta1
    Iota_2 = obj.freq(2,k).*time(:,1,k);                                                                                % Get hyper-time coordinates theta2
     
    Iota_3(:,:,k) = mod(repmat(Iota_1,1,n_char)+obj.phi(1,:),2*pi);                                                     % Add the phase difference of the starting points of the charcteristics and remap values into square [0,2+pi]
    Iota_4(:,:,k) = mod(repmat(Iota_2,1,n_char)+obj.phi(2,:),2*pi);                                                     % Add the phase difference of the starting points of the charcteristics and remap values into square [0,2+pi]
    

    Iota_11(:,k) = reshape(Iota_3(:,:,k),[n_char*reso,1]);                                                              % Reshape data to column vectors
    Iota_22(:,k) = reshape(Iota_4(:,:,k),[n_char*reso,1]);                                                              % Reshape data to column vectors
    Iota_33(:,:,k) = reshape(permute(s0(:,:,:,k),[1,3,2,4]),[n_char*reso,dim]);                                         % Reshape data to column vectors
    
    %Interpolation
    for jk = 1:obj.n
        I3 = Iota_33(:,jk,k);                                                                                           % Get data for interpolation
        F = scatteredInterpolant(Iota_11(:,k),Iota_22(:,k),I3,'natural'); %%                                            % Interpolate scattered data from time integration
        S(:,:,jk,counter) = F(xq,yq);
    end
    clear Iota_1 Iota_2 Iota_3 Iota_4                                                                                   % Clear values
    hypertime  = [reshape(repmat(tmp.',[1,options.resolution]),[1,options.resolution^2]);                               % Set the values for the hypertime coodinates
                  reshape(repmat(tmp  ,[options.resolution,1]),[1,options.resolution^2])];
end                                                                                                                     
s_hypertime = permute(S,[2,1,3,4]);                                                                                     % Define return value hypertime
hypertime = repmat(reshape(hypertime.',[options.resolution,options.resolution,2]),[1,1,1,n_idx]);                       % Define return values of the manifolds
mu = obj.mu(1,index);

% Interpolate boundaries of other torus coordinate
counter = 0;
for k = index
    counter = counter + 1;
    s_hypertime(1,:,:,counter) = 0.5.*(s_hypertime(2,:,:,counter)+s_hypertime(end-1,:,:,counter));                      % Interpolate other boundary which is not smoth due to shooting never set a solution point on the boundary
    s_hypertime(end,:,:,counter) = s_hypertime(1,:,:,counter);                                                          % Set values of opposite boundary equal, due to periodic boundary conditions

    s_hypertime(:,1,:,counter) = 0.5.*(s_hypertime(:,2,:,counter)+s_hypertime(:,end-1,:,counter));                      % Interpolate other boundary which is not smoth due to shooting never set a solution point on the boundary
    s_hypertime(:,end,:,counter) = s_hypertime(:,1,:,counter);                                                          % Set values of opposite boundary equal, due to periodic boundary conditions
end

end
