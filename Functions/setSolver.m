% Method defines ode solver for Shooting algorithms
% prerequisite: obj class has a property called solver_function 
%
%@obj:      Object of a class using shooting algorithm
%@solver:   string defining the solver to be set

function obj = setSolver(obj,solver)

if strcmpi(solver,'ode45')
    obj.solver_function = @(Fcn,T0,z0,options)ode45(Fcn,T0,z0,options);
elseif strcmpi(solver,'ode78')
    obj.solver_function = @(Fcn,T0,z0,options)ode78(Fcn,T0,z0,options);
elseif strcmpi(solver,'ode89')
    obj.solver_function = @(Fcn,T0,z0,options)ode89(Fcn,T0,z0,options);
elseif strcmpi(solver,'ode23')
    obj.solver_function = @(Fcn,T0,z0,options)ode23(Fcn,T0,z0,options);
elseif strcmpi(solver,'ode113')
    obj.solver_function = @(Fcn,T0,z0,options)ode113(Fcn,T0,z0,options);   
elseif strcmpi(solver,'ode15s')
    obj.solver_function = @(Fcn,T0,z0,options)ode15s(Fcn,T0,z0,options);
elseif strcmpi(solver,'ode23s')
    obj.solver_function = @(Fcn,T0,z0,options)ode23s(Fcn,T0,z0,options);
elseif strcmpi(solver,'ode23t')
    obj.solver_function = @(Fcn,T0,z0,options)ode23t(Fcn,T0,z0,options);    
elseif strcmpi(solver,'ode23tb')
    obj.solver_function = @(Fcn,T0,z0,options)ode23tb(Fcn,T0,z0,options);          
end

end