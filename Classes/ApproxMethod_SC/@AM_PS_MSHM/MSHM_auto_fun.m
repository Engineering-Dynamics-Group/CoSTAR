%Function for generating the residuum by means of a periodic shooting method for autonomous function
%
%@obj: ApproxMethod subclass AM_PS_SHM object
%@y:   solution vector for the continuation (contains continuation parameter)
%@DYN  DynamicalSystem class object
%
%@res: resdiual vector for the newton-type corrector method in Continuation class


function res = MSHM_auto_fun(obj,y,DYN)

    % This preallocation of the variable makes the code for some reasons
    % way faster... don't know why...

    Fcn     = DYN.rhs;
    s0      = obj.iv(1:(end-1)); 
    x       = y(1:(end-1));
    mu      = y(end);
    s       = x(1:(end-1),:);
    
    %Evaluate the active parameter (for some reason preallocating these variable is way faster
    % than using them directly)
    param = DYN.param;
    param{DYN.act_param} = mu;

    [~,temp] = obj.solver_function(@(t,y)FCNwrapper(t,y,@(tau,z)Fcn(tau,z,param)),[0,2*pi],x,obj.odeOpts);
    f = Fcn(0,s0,param);

    res = temp(end,1:(end-1)).'-s;
    res(end+1,:) = f.'*(s-s0);              %This is the Poincare condition
    
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






