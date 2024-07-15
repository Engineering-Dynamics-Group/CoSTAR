%Function for generating the residuum by means of a periodic shooting method for non-autonomous function
%
%@obj: ApproxMethod subclass AM_PS_SHM object
%@y:   solution vector for the continuation (contains continuation parameter)
%@DYN  DynamicalSystem class object
%
%@res: resdiual vector for the newton-type corrector method in Continuation class


function res = MSHM_fun(obj,y,DYN)

    Fcn = DYN.rhs;
    mu = y(end);
    n_shoot = obj.n_shoot;                                                  % Number of shooting points
    dim = DYN.dim;
    T = 2*pi./DYN.non_auto_freq(mu);

    dT = T/n_shoot;                                                         % Time for each shooting operation
    Zend = zeros(dim,n_shoot);                                              % Initialize array for end/intermediate points of shooting
    y_perm = zeros(dim*n_shoot,1);                                          % Initialize permuted vector y_perm

    T0 = [0:dT:(n_shoot-1)*dT;dT:dT:n_shoot*dT].';                          % Define time intervals according to number of shooting points
    z0 = reshape(y(1:end-1),[dim,n_shoot]);                                 % reshape y to state-space dimension x number of shooting points 

    Z = cell(n_shoot,1);                                                    % Initialize cell-array for solutions of shooting (necessary 

    %Evaluate the active parameter (for some reason preallocating these variable is way faster
    % than using them directly)
 
    param = DYN.param;
    param{DYN.act_param} = mu;

    for k=1:n_shoot
        [~,Z{k,1}] = obj.solver_function(@(t,z)Fcn(t,z,param),T0(k,:),z0(:,k),obj.odeOpts); 
        Zend(:,k) = Z{k,1}(end,:);
    end
        
    yend = reshape(Zend,[dim*n_shoot,1]);                                   % Reshape Zend to dim*n_shoot x 1 array
    y_perm(1:end-dim,1) = y(dim+1:end,1);                                   % Shift initial vector y to allow direct subtraction
    y_perm(end-dim+1:end,1) = y(1:dim,1);                                   % Shift initial vector y to allow direct subtraction

    res = yend-y_perm;                                                      % Calculate residual

end


