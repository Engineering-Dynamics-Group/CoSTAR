% Method getIV generates the initial value for the defined ApproxMethod method
% from the provided initial condition
function [obj] = getIV(obj,DYN)

if ~isempty(DYN.auto_freq)                                                  %system is (partly) autonomous
    if(obj.n_shoot>1)
        Fcn = DYN.rhs;

        y = DYN.opt_init.ic(:);
        mu = y(end);
        n_shoot = obj.n_shoot;                                                  % Number of shooting points
        dim = DYN.dim;
        T = 2*pi./DYN.auto_freq;

        T0 = linspace(0,T,n_shoot+1);
        param = DYN.param;
        param{DYN.act_param} = mu;
        
        [~,Z] = obj.solver_function(@(t,z)Fcn(t,z,param),T0,y,obj.odeOpts); 

        obj.iv = [reshape(Z(1:end-1,:).',[dim*n_shoot,1]);DYN.auto_freq(:)];
    
    else 
        obj.iv = [DYN.opt_init.ic(:);DYN.auto_freq(:)];                         %Append the initial_condition vector, if it is an autonomous system.
    end
    
else
    
    if(obj.n_shoot>1)
        Fcn = DYN.rhs;

        y = DYN.opt_init.ic(:);
        mu = DYN.non_auto_freq(DYN.param{DYN.act_param});
        n_shoot = obj.n_shoot;                                                  % Number of shooting points
        dim = DYN.dim;
        T = 2*pi./DYN.non_auto_freq(mu);

        T0 = linspace(0,T,n_shoot+1);
        param = DYN.param;
        param{DYN.act_param} = mu;
        
        [~,Z] = obj.solver_function(@(t,z)Fcn(t,z,param),T0,y,obj.odeOpts); 

        obj.iv = reshape(Z(1:end-1,:).',[dim*n_shoot,1]);
    else   
        obj.iv = DYN.opt_init.ic(:);
    end

end

end