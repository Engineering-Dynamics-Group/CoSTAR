%Function for generating the residuum by means of a periodic shooting method for autonomous function
%
%@obj: ApproxMethod subclass AM_PS_SHM object
%@y:   solution vector for the continuation (contains continuation parameter)
%@DYN  DynamicalSystem class object
%
%@res: residual vector for the newton-type corrector method in Continuation class


function res = SHM_auto_fun(obj,y,DYN)

    % This preallocation of the variable makes the code for some reasons
    % way faster... don't know why...
    Fcn     = DYN.rhs;
    n_shoot = obj.n_shoot;
    dim     = DYN.dim;
    s0      = obj.iv(1:(end-1)); 
    x       = y(1:(end-1));
    mu      = y(end);
    s       = x(1:(end-1),:);
    omega   = y(end-1,1);

    dT = 2*pi/n_shoot;                                                      % Time for each shooting operation
    Zend = zeros(dim,n_shoot);                                              % Initialize array for end/intermediate points of shooting
    y_perm = zeros(dim*n_shoot,1);                                          % Initialize permuted vector y_perm
    Z = cell(n_shoot,1);                                                    % Initialize cell-array for solutions of shooting (necessary 
    
    param = DYN.param;
    param{DYN.act_param} = mu;                                              % Evaluate the active parameter
    
    T0 = [0:dT:(n_shoot-1)*dT;dT:dT:n_shoot*dT].';                          % Define time intervals according to number of shooting points  
    z0 = reshape(s,[dim,n_shoot]);                                          % reshape y to state-space dimension x number of shooting points 

  
    f_temp = zeros(dim,n_shoot);
    for k=1:n_shoot
        [~,Z{k,1}] = obj.solver_function(@(t,y)FCNwrapper(t,y,@(tau,z)Fcn(tau,z,param)),T0(k,:),[z0(:,k);omega],obj.odeOpts); 
        Zend(:,k) = Z{k,1}(end,1:dim);
        f_temp(:,k) = Fcn(T0(k,1),s0((k-1)*dim+1:k*dim,1),param);
    end
       
    f = reshape(f_temp,[dim*n_shoot,1]);
    yend = reshape(Zend,[dim*n_shoot,1]);                                   % Reshape Zend to dim*n_shoot x 1 array
    y_perm(1:end-dim,1) = s(dim+1:end,1);                                   % Shift initial vector y to allow direct subtraction
    y_perm(end-dim+1:end,1) = s(1:dim,1);                                   % Shift initial vector y to allow direct subtraction

    res = yend-y_perm; 
    res(end+1,:) = f.'*(s-s0);                                              % This is the Poincare condition

end

% For the autonomous case a function wrapper is needed, which adds the equation T' = 0, 
% which needs to be added to the function and multiplies the right hand side f of the original ODE with T/2*pi
function dzdt = FCNwrapper(t,y,Fcn)
    z = y(1:(end-1),:);
    base_frqn = y(end,:);       
    dzdt = 1./base_frqn.*Fcn(t,z);
    dzdt(end+1,:) = 0;
end
