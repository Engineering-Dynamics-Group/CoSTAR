%Function for generating the residuum by means of a periodic shooting method for non-autonomous function for the calculation of Floquet multipliers 
%
%@obj: Stability subclass ST_PS_SHM object
%@x:   solution vector for the continuation (contains NO continuation parameter)
%@mu:  Continuation parameter
%@DYN  DynamicalSystem class object
%
%@res: resdiual vector for the newton-type corrector method in Stability class

function res = SHM_single_fun(obj,x,mu,DYN)

    Fcn = DYN.rhs;
    T = 2.*pi/(DYN.non_auto_freq(mu));

    %Evaluate the active parameter (for some reason preallocating these variable is way faster
    % than using them directly)
    param = DYN.param;
    param{DYN.act_param} = mu;
    
    [~,temp] = obj.solver_function(@(t,z)Fcn(t,z,param),[0,T],x,obj.odeOpts);
    res = temp(end,:).'-x;

end
