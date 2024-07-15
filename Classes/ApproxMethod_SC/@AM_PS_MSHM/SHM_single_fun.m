%Function for generating the residuum by means of a periodic shooting method for non-autonomous function
%
%@obj: ApproxMethod subclass AM_PS_SHM object
%@y:   solution vector for the continuation (contains continuation parameter)
%@DYN  DynamicalSystem class object
%
%@res: resdiual vector for the newton-type corrector method in Continuation class


function res = SHM_single_fun(obj,y,DYN)

    Fcn = DYN.rhs;
    x = y(1:(end-1));
    mu = y(end);
    T = 2.*pi/(DYN.non_auto_freq(mu));

    %Evaluate the active parameter (for some reason preallocating these variable is way faster
    % than using them directly)
 
    param = DYN.param;
    param{DYN.act_param} = mu;
    
    [~,temp] = obj.solver_function(@(t,z)Fcn(t,z,param),[0,T],x,obj.odeOpts);
    res = temp(end,:).'-x;

end