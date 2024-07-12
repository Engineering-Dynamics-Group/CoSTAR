% This function is a method of the subclass AM_PS_FDM.
% It creates x_IV, which is the major part of the initial value y_IV = [x_IV; mu_IV], for the nonlinear equation solver (e.g. fsolve) to find the initial solution.
% x_IV consists of the method solution vector s_IV and, in the autonomous case, of the autonomous frequency: x_IV = [s_IV; omega_IV]
% Using finite-differences, s_IV contains the state space vectors z(theta,mu_IV) at all discretised points in hyper-time.
% s_IV is created via a 1-st order Fourier series.
%
%@obj:  ApproxMethod subclass AM_PS_FDM object
%@DYN:  DynamicalSystem class object

function obj = getIV(obj,DYN)    

    % Get required parameters
    dim = DYN.system.dim;                       % Dimension of state space
    n_int = obj.n_int;                          % Number of intervals in [0, 2*pi] used for the finte-difference discretization
    

    % If fdm_sol is provided instead of c0, c1 and s1 ...
    if ~isempty(obj.fdm_sol)                                            
        
        n_int_fdm_sol = numel(obj.fdm_sol)/dim; % Get the number of intervals which were used to calculate fdm_sol

        % If n_int_fdm_sol matches n_int ...
        if n_int_fdm_sol == n_int  

            s_IV = obj.fdm_sol;                 % ... fdm_sol can directly be used as s_IV
       
        % If n_int_fdm_sol does not match n_int: Interpolate fdm_sol and evaluate the interpolated data at the theta-values corresponding to n_int
        else

            theta_interp = 0 : 2*pi/n_int_fdm_sol : 2*pi;                   % theta values corresponding to fdm_sol and n_int_fdm_sol
            Z_interp = [reshape(obj.fdm_sol,dim,n_int_fdm_sol), obj.fdm_sol(1:dim)];   % Reshape fdm_sol to Z = [z(0), ... , z(2*pi-DeltaTheta)] and add z(0) (periodic condition)
            Z_csape = csape(theta_interp, Z_interp, 'periodic');            % Interpolation with periodic end conditions (1st and 2nd derivative). Output of csape is a struct (piecewise polynomial form)

            theta_eval = 0 : 2*pi/n_int : 2*pi;                             % theta values corresponding to s_IV and n_int
            s_eval = reshape(fnval(Z_csape, theta_eval), (n_int+1)*dim, 1); % Evaluate the interpolated data and reshape the output to a vector
            s_IV = s_eval(1:n_int*dim);                                     % The first n_int*dim elements are s_IV (last dim elements belong to theta = 2*pi)

        end
    

    % Otherwise: use c0, c1 and s1
    else                                                                
 
        C0 = obj.c0;    if isempty(C0);  C0 = zeros(dim,1);  end        % 0-th order Fourier coefficient
        C1 = obj.c1;    if isempty(C1);  C1 = zeros(dim,1);  end        % 1-st order cosine Fourier coefficient    
        S1 = obj.s1;    if isempty(S1);  S1 = zeros(dim,1);  end        % 1-st order sine Fourier coefficient    
        
        % Create the discretised hypertime vector. The equation system will be build at theta = [0, DeltaTheta, ..., 2*pi-DeltaTheta], which results in n_int equations
        DeltaTheta = 2*pi / n_int;                                      % Hypertime interval between two consecutive disretised points
        theta = 0 : DeltaTheta : (2*pi-DeltaTheta);                     % Create the discretised hyper-time vector. length(theta) = n_int

        % Create a matrix which stores the state space vectors z(theta) of the "start solution": Z = [z(0), ... , z(2*pi-DeltaTheta)]
        Z_IV = repmat(C0,1,n_int) + C1.*cos(theta) + S1.*sin(theta);
     
        % Build up s0 out of Z0
        s_IV = reshape(Z_IV,n_int*dim,1);       % All state space vectors z(theta) are arranged one below the other in the method solution vector s    
   
    end

        
    % Set initial value
    if DYN.n_auto == 0                              
        obj.iv = s_IV;                          % Non-autonomous system: x_IV = s_IV
    elseif DYN.n_auto == 1                      % Autonomous system (periodic solution: n_auto has to be 1)
         obj.iv = [s_IV; DYN.auto_freq];        % Autonomous system. x_IV = [s_IV; omega_IV]
    end

end