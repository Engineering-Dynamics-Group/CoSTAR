% This function is a method of subclass AM_PS_FDM.
% It is needed for comuting the stability via the shooting method and returns the state space vector z(theta=0).
% This vector is used as an initial condition by the shooting solver.
%
%@obj:  ApproxMethod subclass AM_PS_FDM object
%@y:    Solution vector of the continuation in method space (method solution vector, autonomous frequencies, continuaton parameter)
%@DYN:  DynamicalSystem class object
%@IC:   Initial condition vector in state space

function IC = getIC(obj,y,DYN,n_shoot)                                            

    dim = DYN.dim;                                                          % Dimension of state space 
    Omega = DYN.non_auto_freq(y(end,1));
    n_int = obj.n_int;
    T_p = 2.*pi/Omega;
    T = linspace(0,T_p,n_int+1);
    Z = reshape(y(1:end-1,1),[dim,n_int]);

    spl = csape(T,[Z,Z(:,1)],'periodic');                                   % Spline interpolation of FD solution with periodic boundary conditions

    Tint = linspace(0,T_p,n_shoot+1);                                       % Define time interval of shooting points
    Z_int = fnval(spl,Tint);                                                % Evaluate interpolation

    IC = reshape(Z_int(:,1:end-1),[dim*n_shoot,1]);                         % The initial condition vector is the state space vector at the shooting points
 
end