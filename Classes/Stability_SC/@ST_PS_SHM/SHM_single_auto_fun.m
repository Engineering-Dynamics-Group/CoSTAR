%Function for generating the residuum by means of a periodic shooting method for non-autonomous function for the calculation of Floquet multipliers 
%
%@obj: Stability subclass ST_PS_SHM object
%@x:   solution vector for the continuation (contains no continuation parameter)
%@mu:  S
%@DYN  DynamicalSystem class object
%
%@res: resdiual vector for the newton-type corrector method in Stability class


function res = SHM_single_auto_fun(obj,x,mu,x0,DYN)

    % This preallocation of the variable makes the code for some reasons
    % way faster... don't know why...

    Fcn     = DYN.rhs;              %ODE
    s0      = x0(1:(end-1));        %initial value needed for the poincare condition "-2" for autonomous frequency AND continuation parameter.
    s       = x(1:(end-1),:);       %State space vector without autonomous frequency
    
    %Evaluate the active parameter (for some reason preallocating these variable is way faster
    % than using them directly)
    param = DYN.param;
    param{DYN.act_param} = mu;

    [~,temp] = obj.solver_function(@(t,x)FCNwrapper(t,x,@(tau,z)Fcn(tau,z,param)),[0,2*pi],x,obj.odeOpts);
    f = Fcn(0,s0,param);

    res = temp(end,1:(end-1)).'-s;
    res(end+1,:) = f.'*(s-s0);              %This is the Poincare condition
    
end

%For the autonomous case a function wrapper is needed, which adds the
%equation T' = 0, which needs to be added to the function and multiplies the right hand side f of
%the original ODE with T/2*pi
function dzdt = FCNwrapper(t,x,Fcn)
    z = x(1:(end-1),:);
    base_frqn = x(end,:);       
    dzdt = 1./base_frqn.*Fcn(t,z);
    dzdt(end+1,:) = 0;
end






