% This function is a method of the subclass AM_QPS_FDM.
% It creates x_IV, which is the major part of the start vector y0 = [x_IV; mu_IV], for the nonlinear equation solver (e.g. fsolve) to find the initial solution.
% x_IV consists of the method solution vector s_IV and, in the autonomous case, of the autonomous frequency(s): x_IV = [s_IV; (omega_1_IV;) omega_2_IV]
% Using finite-differences, s_IV contains the state space vectors z(theta_1,theta_2,mu_IV) at all discretised points in hyper-time.
% s_IV is created via a multi-dimensional 1-st order Fourier series.
%
%@obj:  ApproxMethod subclass AM_QPS_FDM object
%@DYN:  DynamicalSystem class object

function obj = getIV(obj,DYN)
    
    % Get required parameters
    dim = DYN.system.dim;                       % Dimension of state space
    n_int_1 = obj.n_int_1;                      % Number of intervals in theta_1 direction needed to build the interval [0, 2*pi]
    n_int_2 = obj.n_int_2;                      % Number of intervals in theta_2 direction needed to build the interval [0, 2*pi]


    % If fdm_sol is provided instead of c0, c1_matrix and s1_matrix ...
    if ~isempty(obj.fdm_sol)                                            
        
        % Get the number of intervals which were used to calculate fdm_sol
        if ~isempty(obj.n_int_1_fdm_sol) && ~isempty(obj.n_int_2_fdm_sol)           % If n_int_1_fdm_sol and n_int_2_fdm_sol are provided: Use them
            n_int_1_fdm_sol = obj.n_int_1_fdm_sol;                                  % The Gatekeeper assured that n_int_1_fdm_sol and n_int_2_fdm_sol match fdm_sol
            n_int_2_fdm_sol = obj.n_int_2_fdm_sol;
        elseif isempty(obj.n_int_1_fdm_sol) && isempty(obj.n_int_2_fdm_sol)         % If neither n_int_1_fdm_sol nor n_int_2_fdm_sol are provided: Use n_int_1 and n_int_2
            n_int_1_fdm_sol = n_int_1;                                              % The Gatekeeper assured that n_int_1 and n_int_2 match fdm_sol
            n_int_2_fdm_sol = n_int_2;
        elseif ~isempty(obj.n_int_1_fdm_sol)                                        % If n_int_1_fdm_sol was provided by user, but n_int_2_fdm_sol was not provided
            n_int_1_fdm_sol = obj.n_int_1_fdm_sol;                                  % Use provided n_int_1_fdm_sol ...
            n_int_2_fdm_sol = numel(obj.fdm_sol)/(obj.n_int_1_fdm_sol*dim);         % ... and calculate n_int_2_fdm_sol
        elseif ~isempty(obj.n_int_2_fdm_sol)                                        % If n_int_2_fdm_sol was provided by user, but n_int_1_fdm_sol was not provided
            n_int_1_fdm_sol = numel(obj.fdm_sol)/(obj.n_int_2_fdm_sol*dim);         % Use provided n_int_2_fdm_sol ...
            n_int_2_fdm_sol = obj.n_int_2_fdm_sol;                                  % ... and calculate n_int_1_fdm_sol
        end

        % If n_int_1_fdm_sol matches n_int_1 and n_int_2_fdm_sol matches n_int_2 ...
        if (n_int_1_fdm_sol == n_int_1) && (n_int_2_fdm_sol == n_int_2)

            s_IV = obj.fdm_sol;                 % ... fdm_sol can directly be used as s_IV

        % If the condition above is not fulfilled: Interpolate fdm_sol and evaluate the interpolated data at the theta-values corresponding to n_int_1 and n_int_2
        else

            theta_1_interp = 0 : 2*pi/n_int_1_fdm_sol : 2*pi;           % Create the theta_1-coordinates of fdm_sol ("x1"-values for the interpolation)
            theta_2_interp = 0 : 2*pi/n_int_2_fdm_sol : 2*pi;           % Create the theta_2-coordinates of fdm_sol ("x2"-values for the interpolation)
            Z_reshape = reshape(obj.fdm_sol,[dim,n_int_1_fdm_sol,n_int_2_fdm_sol]);     % Reshape fdm_sol (this structure is needed by csape)
            Z_interp = zeros(dim,n_int_1_fdm_sol+1,n_int_2_fdm_sol+1);                  % Create a new array which will include the state vectors at (theta_1=2*pi, theta_2) and (theta_1, theta_2=2*pi)
            Z_interp(:,1:end-1,1:end-1) = Z_reshape;                    % Fill with Z_reshape
            Z_interp(:,end,1:end-1) = Z_reshape(:,1,:);                 % Add z(theta_1=2*pi, theta_2) = z(theta_1=0, theta_2) for periodic boundary condition (for all theta_2 and for all dimensions of state space)
            Z_interp(:,1:end-1,end) = Z_reshape(:,:,1);                 % Add z(theta_1, theta_2=2*pi) = z(theta_1, theta_2=0) for periodic boundary condition (for all theta_1 and for all dimensions of state space)
            Z_interp(:,end,end) = Z_reshape(:,1,1);                     % Add z(theta_1=2*pi, theta_2=2*pi) = z(theta_1=0, theta_2=0) for periodic boundary condition (for all dimensions of state space)
            Z_csape = csape({theta_1_interp,theta_2_interp},Z_interp,{'periodic','periodic'});  % Interpolation with periodic end conditions (1st and 2nd derivative). Output of csape is a struct (piecewise polynomial form)
            
            theta_1_eval = 0 : 2*pi/n_int_1 : 2*pi*(1-1/n_int_1);       % theta values corresponding to s_IV and n_int_1
            theta_2_eval = 0 : 2*pi/n_int_2 : 2*pi*(1-1/n_int_2);       % theta values corresponding to s_IV and n_int_2
            Z_eval = fnval(Z_csape,{theta_1_eval,theta_2_eval});        % Evaluate the interpolation at the evaluation coordinates {theta_eval_1,theta_eval_2}               
            s_IV = reshape(Z_eval,n_int_1*n_int_2*dim,1);               % Reshape the array

        end
    

    % Otherwise: use c0, c1_matrix and s1_matrix
    else 

        C0 = obj.c0;              if isempty(C0);          C0 = zeros(dim,1);                               end     % 0-th order Fourier coefficient    
        C1_mat = obj.c1_matrix;   if size(C1_mat,2) < 3;   C1_mat = [C1_mat, zeros(dim,3-size(C1_mat,2))];  end     % 1-st order cosine Fourier coefficients    
        S1_mat = obj.s1_matrix;   if size(S1_mat,2) < 3;   S1_mat = [S1_mat, zeros(dim,3-size(S1_mat,2))];  end     % 1-st order sine Fourier coefficients 
  

        % Create the discretised hyper-time vectors. The equation system will be build at theta = [0, DeltaTheta, ..., 2*pi-DeltaTheta], which results in n_int equations
        DeltaTheta_1 = 2*pi / n_int_1;                                      % Hyper-time interval between two consecutive disretised points in theta_1 direction
        theta_1 = 0 : DeltaTheta_1 : (2*pi-DeltaTheta_1);                   % Create the discretised hyper-time vector theta_1  ->  length(theta_1) = n_int_1
        Theta_1 = repmat(theta_1,1,n_int_2);                                % Theta_1 = [theta_1, ..., theta_1] <n_int_2-times> ->  length(Theta_1) = n_int_1 * n_int_2
    
        DeltaTheta_2 = 2*pi / n_int_2;                                      % Hyper-time interval between two consecutive disretised points in theta_2 direction
        theta_2 = 0 : DeltaTheta_2 : (2*pi-DeltaTheta_2);                   % Create the discretised hyper-time vector theta_2  ->  length(theta_2) = n_int_2
        Theta_2 = reshape(repmat(theta_2,n_int_1,1),1,n_int_1*n_int_2);     % Theta_2 = [theta_2(1) <n_int_1-times>, ... , theta_2(end) <n_int_1-times>]  ->  length(Theta_2) = n_int_1 * n_int_2


        % Create a matrix which stores the state space vectors z(theta_1,theta_2) for the initial value for fsolve
        % The state space vectors are arranged as follows: Z = [z(theta_1,theta_2(1)), z(theta_1,theta_2(2)), ... z(theta_1,theta_2(end))]
        % -> Z = [z(0,0), ... , z(theta_1(end),0), z(0,DeltaTheta_2), ... , z(theta_1(end),DeltaTheta_2), ...... , z(theta_1(end),theta_2(end))]
        % -> theta_1 is gone through at theta_2 = 0 -> theta_2 is iterated to theta_2(2) and theta_1 is gone through again -> theta_2 is iterated and so on
        % This arrangement does not require a for-loop and it is the reason why Theta_1 and Theta_2 need to have their special structure
        Z_IV = repmat(C0,1,n_int_1*n_int_2) + C1_mat(:,1).*cos(Theta_1) + C1_mat(:,2).*cos(Theta_2) + C1_mat(:,3).*cos(Theta_1+Theta_2) ...
                                            + S1_mat(:,1).*sin(Theta_1) + S1_mat(:,2).*sin(Theta_2) + S1_mat(:,3).*sin(Theta_1+Theta_2);
    
        % Plot the 2-dimensional surface of Z0(1,:) in hyper-time domain (theta_1,theta_2) to check if Z0 is alright
        %{
        [T1,T2] = meshgrid(theta_1,theta_2);
        Z1 = reshape(Z0(3,:),n_int_1,n_int_2)';
        urf(T1,T2,Z1)
        xlabel('\theta_1'); ylabel('\theta_2'); zlabel('Z0_1')
        %}
    
        % Build up s_IV out of Z0
        s_IV = reshape(Z_IV,n_int_1*n_int_2*dim,1);     % All state space vectors z(theta_1,theta_2) are arranged one below the other according to the arrangement scheme described above
    
    end

    
    % Set initial value
    if DYN.n_auto == 0                          % Non-autonomous system (0 autonomous frequencies and 2 non-autonomous frequencies)    
        obj.iv = s_IV;                          % x_IV = s_IV
    elseif DYN.n_auto == 1                      % Partly autonomous system (1 autonomous frequency and 1 non-autonomous frequency)
         obj.iv = [s_IV; DYN.auto_freq];        % x_IV = [s_IV; omega_2_IV]
    elseif DYN.n_auto == 2                      % Full autonomous system (2 autonomous frequencies and 0 non-autonomous frequencies)
        omega_1 = DYN.auto_freq(1);             % Get the first autonomous frequency
        omega_2 = DYN.auto_freq(2);             % Get the second autonomous frequency
        obj.iv = [s_IV; omega_1; omega_2];      % x_IV = [s_IV; omega_1_IV; omega_2_IV]. DYN.auto_freq is not directly inserted into x_IV, because it can be a column OR row vector
    end

end