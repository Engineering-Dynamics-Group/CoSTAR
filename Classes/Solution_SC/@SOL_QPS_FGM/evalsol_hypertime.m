% Method of SOL_QPS_FGM: This method calculates the hypertime manifolds of solution, which is 2-dimensional
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

    %% Parameter
    index = options.index;
    n_evals = numel(index);                                             % Number of solutions to evaluate
    dim = DYN.dim;                                                      % Dimension of the state space
    if size(options.resolution,2) == 2
        res_1 = options.resolution(1);                                  % Desired resolution of output s in theta_1-direction if user defined [1x2] array
        res_2 = options.resolution(2);                                  % Desired resolution of output s in theta_2-direction if user defined [1x2] array
    else
        res_1 = options.resolution;                                     % Desired resolution of output s in theta_1-direction if user defined scalar
        res_2 = res_1;                                                  % Desired resolution of output s in theta_2-direction if user defined scalar
    end
   

    %% Initialize
    s_hypertime = zeros(res_1,res_2,dim,n_evals);                       % Array for the solution output
    mu = obj.mu(1,index);                                               % Get the mu-values
    s = obj.s(1,index);                                                 % Get solution vectors at the desired index
    hmatrix = obj.hmatrix(1,index);                                     % Get the hmatrix at the desired index
    n_hh = obj.n_hh(1,index);                                           % Get the number of higher harmonics at the desired index

    % Compute the theta-values for the evaluation
    theta_eval_1 = linspace(0, 2*pi, res_1);                            % Coordinates of the output in theta_1-direction
    theta_eval_2 = linspace(0, 2*pi, res_2);                            % Coordinates of the output in theta_2-direction
    [Theta_2, Theta_1] = meshgrid(theta_eval_2,theta_eval_1);           % Theta_2 is the first output of meshgrid so that hypertime(:,:,2) has the desired structure
    hypertime(:,:,1) = Theta_1;                                         % Theta_1 values as a matrix. First row is the zero vector, second row only consists of values 2*pi/res etc.
    hypertime(:,:,2) = Theta_2;                                         % Theta_2 values as a matrix. First column is the zero vector, second column only consists of values 2*pi/res etc.
    hypertime = repmat(hypertime,[1,1,1,n_evals]);                      % The matrices are the same for each evaluation
    hypertime_eval = [reshape(Theta_1,1,res_1*res_2);                   % [2 x res_1*res_2] matrix of all (theta_1,theta_2) combinations. This is needed for the evaluation
                      reshape(Theta_2,1,res_1*res_2)];


    %% Evaluate each solution
    for k = 1:n_evals

        s_k = s{1,k};                                                   % Get the vector of Fourier coefficients      
        hmatrix_k = hmatrix{1,k};                                       % Get the corresponding hmatrix
        n_hh_k = n_hh(k);                                               % Get the corresponding number of harmonics

        %This could also be done in the real domain...however: this is how we do it in the residuum equation
        FC = [s_k(1:dim,1);s_k((dim+1):(n_hh_k)*dim,1)-1i.*s_k(((n_hh_k)*dim+1):end,1)];        % Assemble complex Fourrier vector    
        FCtemp = reshape(FC,dim,n_hh_k);
 
        % How does this work?: I am evaluating here the arguments in the cosine or sine functions: cos(H1*theta1+H2*theta2). hmatrix'*p_freq_val is the scalar product
        % The definition of p_freq_val (similar to the meshgrid function) ensures, that p_chf is evaluated at every point of the discretized torus         
        chf = exp(1i.*hmatrix_k'*hypertime_eval);
        % Complex harmonic functions:
        %  - real(obj.p_chf(2,:)) gives the cosine with the first higher harmonic defined in obj.hmatrix (first one is 0)
        %  - imag(obj.p_chf(2,:)) correspondingly gives the same for the sine func                        
   
        s_hypertime(:,:,:,k) = reshape(real(FCtemp*chf).',[res_1,res_2,dim]);  

    end
    
end