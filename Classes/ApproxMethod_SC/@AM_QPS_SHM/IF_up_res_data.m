% IF_up_res_data for Quasiperiodic Shooting
% This function is a method of subclass AM_QPS_SHM. 
% This method updates the necessary data of the ApproximationMethod object 
% with data from the Continuation object
%
% @obj:  ApproximationMethod subclass object
% @CON:  Continuation class object
% @DYN:  DynamicalSystem class object
%
% @obj:  ApproximationMethod object
%
function obj = IF_up_res_data(obj,var1,DYN)

if isa(var1,'Continuation')                                                         % If var1 is an object of Continuation
    obj.iv = var1.yp(1:(end-1));                                                    % Set iv to predictor point
    y_old  = var1.p_y0_old{end};                                                    % Get last curve point
elseif isa(var1,'double')                                                           % var1 should be a solution vector (type double) in all other cases
    obj.iv = var1;                                                                  % Set iv to given solution vector y0
    y_old  = var1;                                                                  % Use given solution
end


%% If autonomous frequency present
% Calculate reference solution and evaluate the derivative of the reference
% solution to be able to evaluate the phase condition(s)
if(DYN.n_auto==1)
    dim = DYN.dim;                                                                  % Get dimension of state-space
    n_char = obj.n_char;                                                            % Get number of charecteristics
    reso = obj.reso_phase;                                                          % Get resolution for time-integration
    
    param = DYN.param;                                                              % Set parameter vector
    param{DYN.act_param} = y_old(end,1);                                            % Set active parameter
    Xchar = linspace(0,2*pi*(1-1/n_char),n_char);                                   % Define interval for gradient
    IV = reshape(y_old(1:end-2,1),[dim,n_char]);                                    % Fetch initial conditions
    Omega = [DYN.non_auto_freq(y_old(end,1)),y_old(end-1,1)];                       % Set Omega according to autonomous frequency and bifurcation parameter
    T_char = linspace(0,2*pi/Omega(1,1),reso);                                      % Integration time for characteristics
    F1 = zeros(reso,n_char,dim);                                                    % Initialize matrix for time solution
        
    PHI(1,:) = 1./Omega(1,1).*obj.phi(1,:);                                         % Define phase for characteristics
    PHI(2,:) = 1./Omega(1,2).*obj.phi(2,:);                                         
    
    FW = @(t,z)obj.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                % Function wrapper to integrate over all characteristics at once
    [~,V] = obj.solver_function(FW,T_char,IV,obj.odeOpts);                          % Time integration of all characteristics
    W = permute(reshape(V,[reso,dim,n_char]),[1,3,2]);                              % Reshape solution
  
    for j = 1:dim                                                                   % Calculate derivative of reference solution for phase condition
        [~,F1(:,:,j)] = gradient(W(:,:,j),obj.phi(2,2),Xchar(1,2));
    end

    obj.Y_old{1,1} = Xchar;                                                         % Save integration interval
    obj.Y_old{1,2} = W;                                                             % Save reference solution
    obj.Y_old{1,3} = F1;                                                            % Save derivative of reference solution with respect to theta2
    
elseif(DYN.n_auto==2)
    dim = DYN.dim;                                                                  % Get dimension of state-space
    n_char = obj.n_char;                                                            % Get number of charecteristics                         
    reso = obj.reso_phase;                                                          % Get resolution for time-integration
    
    param = DYN.param;                                                              % Set parameter vector
    param{DYN.act_param} = y_old(end,1);                                            % Set active parameter
    Xchar = linspace(0,2*pi*(1-1/n_char),n_char);                                   % Define interval for gradient
    IV = reshape(y_old(1:end-3,1),[dim,n_char]);                                    % Fetch initial conditions
    Omega = [y_old(end-2,1),y_old(end-1,1)];                                        % Set Omega according to autonomous frequency and bifurcation parameter
    T_char = linspace(0,2*pi/Omega(1,1),reso);                                      % Integration time for characteristics
    F1 = zeros(reso,n_char,dim);                                                    % Initialize matrix for time solution
    F2 = zeros(reso,n_char,dim);                                                    % Initialize matrix for time solution

    PHI(1,:) = 1./Omega(1,1).*obj.phi(1,:);                                         % Define phase for characteristics
    PHI(2,:) = 1./Omega(1,2).*obj.phi(2,:);                                         
    
    FW = @(t,z)obj.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                % Function wrapper to integrate over all characteristics at once
    [~,V] = obj.solver_function(FW,T_char,IV,obj.odeOpts);                          % Time integration of all characteristics
    W = permute(reshape(V,[reso,dim,n_char]),[1,3,2]);                              % Reshape solution
  
    for j = 1:dim                                                                   % Calculate derivatives of reference solution for phase condition
        [F2(:,:,j),F1(:,:,j)] = gradient(W(:,:,j),obj.phi(2,2),Xchar(1,2));
    end
    obj.Y_old{1,1} = Xchar;                                                         % Save integration interval
    obj.Y_old{1,2} = W;                                                             % Save reference solution
    obj.Y_old{1,3} = F1;                                                            % Save derivative of reference solution with respect to theta1
    obj.Y_old{1,4} = F2;                                                            % Save derivative of reference solution with respect to theta2
end

end