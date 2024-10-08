%Function for generating the residuum by means of a periodic shooting method for autonomous function
%
%@obj: ApproxMethod subclass AM_PS_SHM object
%@y:   solution vector for the continuation (contains continuation parameter)
%@DYN  DynamicalSystem class object
%
%@res: resdiual vector for the newton-type corrector method in Continuation class


function res = MSHM_auto_fun(obj,x,x0,mu,DYN)

    Fcn     = DYN.rhs;
    n_shoot = obj.n_shoot;                                                  % Number of shooting points
    dim = DYN.dim;

    Zend = zeros(dim+1,n_shoot);                                            % Initialize array for end/intermediate points of shooting
    Z = cell(n_shoot,1);
    y_perm = zeros(dim*n_shoot,1);                                          % Initialize permuted vector y_perm

    s0 = x0(1:end-1,1);
    s  = x(1:(end-1),:);
    dT = 2*pi/n_shoot;
    T0 = [0:dT:(n_shoot-1)*dT;dT:dT:n_shoot*dT].';                          % Define time intervals according to number of shooting points
    z0 = [reshape(s,[dim,n_shoot]);repmat(x(end,1),[1,n_shoot])];           % reshape s to state-space dimension x number of shooting points and add autonomous frequency

    param = DYN.param;
    param{DYN.act_param} = mu;

    for k=1:n_shoot
        [~,Z{k,1}] = obj.solver_function(@(t,y)FCNwrapper(t,y,@(tau,z)Fcn(tau,z,param)),T0(k,:),z0(:,k),obj.odeOpts); 
        Zend(:,k) = Z{k,1}(end,:);
    end
    f = Fcn(0,s0,param);
    
    yend = reshape(Zend(1:end-1,:),[dim*n_shoot,1]);                        % Reshape Zend to dim*n_shoot x 1 array
    y_perm(1:end-dim,1) = s(dim+1:end,1);                                   % Shift initial vector y to allow direct subtraction
    y_perm(end-dim+1:end,1) = s(1:dim,1);                                   % Shift initial vector y to allow direct subtraction

    res = yend-y_perm;                                                      % Calculate residual
    res(end+1,:) = f.'*(s(1:dim,1)-s0(1:dim,1));                            % This is the Poincare condition
    
end

%For the autonomous case a function wrapper is needed, which adds the
%equation T' = 0, which needs to be added to the function and multiplies the right hand side f of
%the original ODE with T/2*pi
function dzdt = FCNwrapper(t,y,Fcn)
    z = y(1:(end-1),:);
    base_frqn = y(end,:);       
    dzdt = 1./base_frqn.*Fcn(t,z);
    dzdt(end+1,:) = 0;
end






