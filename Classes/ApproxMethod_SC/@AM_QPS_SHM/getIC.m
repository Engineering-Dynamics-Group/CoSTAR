% This function is a method of subclass AM_QPS_SHM.
% It is needed for computing the stability when n_char ~= n_char_st.
% If the above statement is true, the initial state vectors z(0,theta_2) are interpolated
%
% @obj: ApproxMethod subclass AM_QPS_SHM object
% @y:   Solution vector of the continuation (method solution vector, autonomous frequencies, continuaton parameter)
% @DYN: DynamicalSystem class object
% @IC:  Initial condition vector in state space

function IC = getIC(obj,y,DYN,n_char_st)                                            

    % Parameters
    dim = DYN.dim;                                                      % Dimension of state space 
    n_auto = DYN.n_auto;                                                % Number of autonomous frequencies
    s = y(1:end-1-n_auto);                                              % Method solution vector
    n_char = obj.n_char;                                                % Number of characteristics of quasi-periodic shooting method

    % Get the initial state vectors z(0,theta_2)
    Z0 = reshape(s,dim,n_char);                                         % Reshape the computed initial states to a [dim x n_char] matrix for interpolation
    theta_2 = linspace(0,2*pi*(1-1/n_char),n_char);                     % theta_2 values of the computed solution
    
    % Interpolation
    Z_interp = csape([theta_2,2*pi],[Z0,Z0(:,1)],'periodic');           % Interpolate the initial states (do not forget to append by the first element to have periodic boundary)
    
    % Evaluation
    theta_2_eval = linspace(0,2*pi*(1-1/n_char_st),n_char_st);          % theta_2 values for the new number of characteristics n_char_st
    Z_eval = fnval(Z_interp, theta_2_eval);                             % Evaluate interpolation at the desired theta_2 values
    
    % Initial states for n_char_st
    IC = reshape(Z_eval,dim*n_char_st,1);                               % Reshape the new initial states to a vector

end