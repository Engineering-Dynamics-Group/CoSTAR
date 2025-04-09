% Method getIV for Quasiperiodic Shooting
% This function is a method of subclass AM_QPS_SHM.
% This method generates the initial value for qps shooting method from the
% provided initial condition in state-space if the solution is stable
%
% @obj:  ApproximationMethod subclass object
% @DYN:  DynamicalSystem class object
%
% @obj:  ApproximationMethod object
%
function obj = getIV(obj,DYN)

%% Initialization
obj.mu0 = DYN.param{DYN.act_param};                                                             % Set obj.mu0 to active parameter in DYN.param
mu0 = obj.mu0;                                                                                  % Get mu0
reso = obj.reso_phase;                                                                          % Get resolution of solution
dim = DYN.dim;                                                                                  % Get dimension of state-space
n_char = obj.n_char;                                                                            % Get number of characteristics
F1 = zeros(reso,n_char,dim);                                                                    % Initialize F1 and F2 for derivative of reference solution
F2 = zeros(reso,n_char,dim);


%% If user supplied initial value
if isfield(DYN.opt_init,'iv')

    param = DYN.param;                                                                          % Get param vector
    param{DYN.act_param} = mu0;                                                                 % Set param vector to mu0
    iv = reshape(DYN.opt_init.iv,[],1);                                                         % Supplied initial value. Make it a column vector if it is a row vector
    n_char_iv =  numel(iv)/dim;                                                                 % n_char of iv
    if n_char_iv ~= n_char                                                                      % If iv was calculated using a different number of characteristics (~= n_char)
        phi_0_iv = linspace(0,2*pi*(1-1/n_char_iv),n_char_iv);                                  % Define the spacing of the characteristics of iv
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        iv_interp = csape([phi_0_iv,2*pi],[reshape(iv,dim,n_char_iv),iv(1:dim)],'periodic');    % Interpolate initial values
        IV = reshape(fnval(iv_interp,phi_0),dim*n_char,1);                                      % Evaluate interpolated iv at phi_0 and reshape output to a vector
    else                                                                                        % If iv was calculated using n_char
        IV = iv;                                                                                % Set iv to IV
    end

    if DYN.n_auto==0
        Omega = DYN.non_auto_freq(mu0);
        if(Omega(1,1)>Omega(1,2)); s1=1; s2=0; else; s1=0; s2=1; end                            % Check which "direction" for time-integration is best
        T = s1.*2*pi/Omega(1,1) + s2.*2*pi/Omega(1,2);                                          % Define period
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        PHI(1,:) = s1.*obj.phi(1,1)*ones(1,length(phi_0)) + s2.*(obj.phi(1,1) + phi_0);         % Generate the values for the phase shift for either excitation
        PHI(2,:) = s2.*obj.phi(2,1)*ones(1,length(phi_0)) + s1.*(obj.phi(2,1) + phi_0);
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    
    elseif DYN.n_auto==1
        Omega(1,1) = DYN.non_auto_freq(mu0);                                                    % Set Omega1 to non-autonomous frequency
        Omega(1,2) = DYN.auto_freq;                                                             % Set Omega2 to provided autonomous frequency
        T = 2*pi/Omega(1,1);                                                                    % Define Integration period
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    
    elseif DYN.n_auto==2
        Omega(1,1) = DYN.auto_freq(1);                                                          % Set Omega1 to first provided autonomous frequency
        Omega(1,2) = DYN.auto_freq(2);                                                          % Set Omega2 to second provided autonomous frequency
        T = 2*pi/Omega(1,1);                                                                    % Define Integration period
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    
    end

    Xchar = linspace(0,2*pi*(1-1/n_char),n_char);
    obj.Ik = [0,T];

    
%% If user did not supply initial value
else

    param = DYN.param;                                                                          % Get param vector
    param{DYN.act_param} = mu0;                                                                 % Set param vector to mu0

    if DYN.n_auto==0
        Omega = DYN.non_auto_freq(mu0);
        if(Omega(1,1)>Omega(1,2)); s1=1; s2=0; else; s1=0; s2=1; end                            % Check which "direction" for time-integration is best
        T = s1.*2*pi/Omega(1,1) + s2.*2*pi/Omega(1,2);                                          % Define period
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        PHI(1,:) = s1.*obj.phi(1,1)*ones(1,length(phi_0)) + s2.*(obj.phi(1,1) + phi_0);         % Generate the values for the phase shift for either excitation
        PHI(2,:) = s2.*obj.phi(2,1)*ones(1,length(phi_0)) + s1.*(obj.phi(2,1) + phi_0);
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    
    elseif DYN.n_auto==1
        Omega(1,1) = DYN.non_auto_freq(mu0);                                                    % Set Omega1 to non-autonomous frequency
        Omega(1,2) = DYN.auto_freq;                                                             % Set Omega2 to provided autonomous frequency
        T = 2*pi/Omega(1,1);                                                                    % Define Integration period
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    
    elseif DYN.n_auto==2
        Omega(1,1) = DYN.auto_freq(1);                                                          % Set Omega1 to first provided autonomous frequency
        Omega(1,2) = DYN.auto_freq(2);                                                          % Set Omega2 to second provided autonomous frequency
        T = 2*pi/Omega(1,1);                                                                    % Define Integration period
        phi_0 = linspace(0,2*pi*(1-1/n_char),n_char);                                           % Define the spacing of the characteristics
        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing

    end

    Xchar = linspace(0,2*pi*(1-1/n_char),n_char);
    obj.Ik = [0,T];


    C0 = obj.c0;              if isempty(C0);          C0 = zeros(dim,1);                               end     % 0-th order Fourier coefficient
    C1_mat = obj.c1_matrix;   if size(C1_mat,2) < 3;   C1_mat = [C1_mat, zeros(dim,3-size(C1_mat,2))];  end     % 1-st order cosine Fourier coefficients
    S1_mat = obj.s1_matrix;   if size(S1_mat,2) < 3;   S1_mat = [S1_mat, zeros(dim,3-size(S1_mat,2))];  end     % 1-st order sine Fourier coefficients

    % Calculate the initial values from a Fourier series
    theta_1 = zeros(1,n_char);                                                                  % We only need the values for theta_1 = 0 
    theta_2 = linspace(0,2*pi*(1-1/n_char),n_char);                                             % Use n_char points in theta_2 - direction

    % Create a matrix which stores the state space vectors z(theta_1=0,theta_2) for the initial value for fsolve
    % The state space vectors are arranged as follows: Z = [z(0,theta_2(1)), ... , z(0,theta_2(end))]
    Z_IV = repmat(C0,1,n_char) + C1_mat(:,1).*cos(theta_1) + C1_mat(:,2).*cos(theta_2) + C1_mat(:,3).*cos(theta_1+theta_2) ...
                               + S1_mat(:,1).*sin(theta_1) + S1_mat(:,2).*sin(theta_2) + S1_mat(:,3).*sin(theta_1+theta_2);
    
    IV = reshape(Z_IV,[dim*n_char,1]);                                                          % Reshape the matrix to a vector

end



%% Save solution
if DYN.n_auto==0

    obj.iv = IV;


elseif DYN.n_auto==1

    obj.iv = [IV;Omega(1,2)];                                                                   % Append obj.iv with autonomous frequencies
    T_char = linspace(0,2*pi/Omega(1,1),reso);                                                  % Integration time for characteristics

    PHI_s(1,:) = 1./Omega(1,1).*obj.phi(1,:);                                                   % Define spacing for characteristics
    PHI_s(2,:) = 1./Omega(1,2).*obj.phi(2,:);

    % Get reference solution for phase-condition
    FW = @(t,z)obj.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI_s);                          % Function wrapper to integrate all characteristics at once
    [~,V] = obj.solver_function(FW,T_char,IV,obj.odeOpts);                                      % Integrate along characteristics to be able to calculate gradient
    W = permute(reshape(V,[reso,dim,n_char]),[1,3,2]);                                          % Reshape solution
    
    for j = 1:DYN.dim
        [~,F1(:,:,j)] = gradient(W(:,:,j),PHI(2,2),Xchar(1,2));                                 % Calculate gradient
    end
    obj.Y_old{1,1} = Xchar;                                                                     % Save integration interval
    obj.Y_old{1,2} = W;                                                                         % Save reference solution
    obj.Y_old{1,3} = F1;                                                                        % Save derivative of reference solution with respect to theta2


elseif DYN.n_auto==2

    obj.iv = [IV;Omega(1,1);Omega(1,2)];                                                        % Append obj.iv with autonomous frequencies
    T_char = linspace(0,2*pi/Omega(1,1),reso);                                                  % Integration time for characteristics

    % Get reference solution for phase-conditions
    PHI_s(1,:) = 1./Omega(1,1).*obj.phi(1,:);                                                   % Define spacing for characteristics
    PHI_s(2,:) = 1./Omega(1,2).*obj.phi(2,:);

    FW = @(t,z)obj.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI_s);                          % Function wrapper to integrate all characteristics at once
    [~,V] = obj.solver_function(FW,T_char,IV,obj.odeOpts);                                      % Integrate along characteristics to be able to calculate gradient
    W = permute(reshape(V,[reso,dim,n_char]),[1,3,2]);                                          % Reshape solution

    for j = 1:DYN.dim
        [F2(:,:,j),F1(:,:,j)] = gradient(W(:,:,j),PHI(2,2),Xchar(1,2));                         % Calculate derivatives
    end
    obj.Y_old{1,1} = Xchar;                                                                     % Save integration interval
    obj.Y_old{1,2} = W;                                                                         % Save reference solution
    obj.Y_old{1,3} = F1;                                                                        % Save derivative of reference solution with respect to theta1
    obj.Y_old{1,4} = F2;                                                                        % Save derivative of reference solution with respect to theta2

end

end