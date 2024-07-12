% This function is a method of subclass AM_PS_FDM.
% It is needed for comuting the stability via the shooting method and returns the state space vector z(theta=0).
% This vector is used as an initial condition by the shooting solver.
%
%@obj:  ApproxMethod subclass AM_PS_FDM object
%@y:    Solution vector of the continuation in method space (method solution vector, autonomous frequencies, continuaton parameter)
%@DYN:  DynamicalSystem class object
%@IC:   Initial condition vector in state space

function IC = getIC(obj,y,DYN)                                            

    dim = DYN.dim;                                  % Dimension of state space 
    
    IC = y(1:dim);                                  % The initial condition vector is the state space vector z(theta=0)
 
end