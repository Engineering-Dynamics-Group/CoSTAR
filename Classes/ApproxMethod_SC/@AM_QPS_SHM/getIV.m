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

if(isfield(DYN.opt_init,'iv'))                                                                  % If user supplies initial value
    param = DYN.param;                                                                          % Get param vector
    param{DYN.act_param} = mu0;                                                                 % Set param vector to mu0
    iv = reshape(DYN.opt_init.iv,[],1);                                                         % Supplied initial value. Make it a column vector if it is a row vector
    n_char_iv =  numel(iv)/dim;                                                                 % n_char of iv
    if n_char_iv ~= n_char                                                                      % If iv was calculated using a different number of characteristics (~= n_char)
        phi_0_iv = linspace(0,2*pi*(1-1/(n_char_iv+1)),n_char_iv);                              % Define the spacing of the characteristics of iv
        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics
        iv_interp = csape([phi_0_iv,2*pi],[reshape(iv,dim,n_char_iv),iv(1:dim)],'periodic');    % Interpolate initial values
        IV = reshape(fnval(iv_interp,phi_0),dim*n_char,1);                                      % Evaluate interpolated iv at phi_0 and reshape output to a vector
    else                                                                                        % If iv was calculated using n_char
        IV = iv;                                                                                % Set iv to IV
    end
    if(DYN.n_auto==0)
        Omega = DYN.non_auto_freq(mu0);
        if(Omega(1,1)>Omega(1,2)); s1=1; s2=0; else; s1=0; s2=1; end                            % Check which "direction" for time-integration is best
        T = s1.*2*pi/Omega(1,1) + s2.*2*pi/Omega(1,2);                                          % Define period
        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics

        PHI(1,:) = s1.*obj.phi(1,1)*ones(1,length(phi_0)) + s2.*(obj.phi(1,1) + phi_0);         % Generate the values for the phase shift for either excitation
        PHI(2,:) = s2.*obj.phi(2,1)*ones(1,length(phi_0)) + s1.*(obj.phi(2,1) + phi_0);
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    elseif(DYN.n_auto==1)
        Omega(1,1) = DYN.non_auto_freq(mu0);                                                    % Set Omega1 to non-autonomous frequency
        Omega(1,2) = DYN.auto_freq;                                                             % Set Omega2 to provided autonomous frequency
        T = 2*pi/Omega(1,1);                                                                    % Define Integration period
        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics

        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    elseif(DYN.n_auto==2)
        Omega(1,1) = DYN.auto_freq(1);                                                          % Set Omega1 to first provided autonomous frequency
        Omega(1,2) = DYN.auto_freq(2);                                                          % Set Omega2 to second provided autonomous frequency
        T = 2*pi/Omega(1,1);                                                                    % Define Integration period
        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics

        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to calculated spacing
    end
    Xchar = linspace(0,2*pi*(1-1/(n_char+1)),n_char);
    obj.Ik = [0,T];
else
    % If no initial value is supplied take supplied initial condition and
    % do a time integration
    if(isfield(DYN.opt_init,'ic'))
        IC = DYN.opt_init.ic;
    else
        IC = zeros(dim,1);
    end
    
    disp('--------------------------')
    disp('Calculating initial value!')
    disp('--------------------------')
    
    if(DYN.n_auto==0)                                                                           % Non-autonomous case
        Omega = DYN.non_auto_freq(DYN.param{DYN.act_param});                                    % Set Omega to given non-autonomous frequencies
        if(Omega(1,1)>Omega(1,2)); s1=1; s2=0; else; s1=0; s2=1; end                            % Define period
        T = s1.*2*pi/Omega(1,1) + s2.*2*pi/Omega(1,2);                                          % Check which "direction" is best for integration

        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics

        PHI(1,:) = s1.*obj.phi(1,1)*ones(1,length(phi_0)) + s2.*(obj.phi(1,1) + phi_0);         % Generate the values for the phase shift for either excitation
        PHI(2,:) = s2.*obj.phi(2,1)*ones(1,length(phi_0)) + s1.*(obj.phi(2,1) + phi_0);
        obj.phi = PHI;                                                                          % Set obj.phi to caculated spacing

    elseif(DYN.n_auto==1)                                                                       % Mixed case
        Omega = DYN.non_auto_freq(DYN.param{DYN.act_param});                                    % Set Omega to non-autonous frequency
        T = 2*pi/Omega;                                                                         % Set integration period
        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics

        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to caculated spacing

    elseif(DYN.n_auto==2)                                                                       % Full-autonomous case
        Omega(1,1) = DYN.auto_freq(1,1);                                                        % Set Omega to first autonous frequency
        T = 2*pi/Omega;                                                                         % Define intergration interval
        phi_0 = linspace(0,2*pi*(1-1/(n_char+1)),n_char);                                       % Define the spacing of the characteristics

        PHI(1,:) = obj.phi(1,1)*ones(1,length(phi_0));                                          % Generate the values for the phase shift for either excitation
        PHI(2,:) = obj.phi(2,1) + phi_0;
        obj.phi = PHI;                                                                          % Set obj.phi to caculated spacing
    end

    %% Compute initial values
    obj.Ik = [0,T];                                                                             % Set integration interval

    %% Initialization   

    %Set frequencies
    if(DYN.n_auto==0)
        Omega = DYN.non_auto_freq(mu0);
    elseif(DYN.n_auto==1)
        Omega(1,1) = DYN.non_auto_freq(mu0);
        Omega(1,2) = DYN.auto_freq;
    elseif(DYN.n_auto==2)
        Omega(1,2) = DYN.auto_freq(1,2);                                                        % Set second autonomous frequency
    end

    param = DYN.param;                                                                          % Set parameter vector
    param{DYN.act_param} = mu0;                                                                 % Set active parameter
    Fcn = @(t,z)DYN.rhs(t+obj.phi(:,1),z,param);                                                % Set right-hand-side

    %% Time integration   
    Tend = obj.tinit;                                                                           % Set end-time for "transient integration"
    Tfin = Tend + obj.deltat;                                                                   % Set end-time for "stationary time integration"
    dt = obj.dt;                                                                                % Set time increment for interval


    T0 = linspace(Tend,Tfin,round((Tfin-Tend)/dt,0));                                           % Set interval for integration on manifold
    [~,Ztemp] = obj.solver_function(Fcn,[0,Tend],IC,obj.odeOpts);                               % Transient time-integration
    [T,Z] = obj.solver_function(Fcn,T0,Ztemp(end,:).',obj.odeOpts);                             % Stationary time-integration on manifold

    %% Isolation of autonomous frequency
    if(DYN.n_auto==1)
        [f,S,psi] = costarFFT(T,Z(:,1));                                                        % Do a fourier-transformation
        ind_u = find(2.*pi.*f>0.95.*Omega(1,2),1);                                              % Find the frequency which is 5% lower than the estimated value for the autonomous frequency
        ind_o = find(2.*pi.*f>1.05.*Omega(1,2),1);                                              % Find the frequency which is 5% larger than the estimated value for the autonomous frequency
        f_reduced = f(ind_u:ind_o,1);                                                           % Make a window in the frequency domain
        S_reduced = S(ind_u:ind_o,1);
        [A,B] = sort(S_reduced,1,'descend');                                                    % Sort the frequencies in desceding order of amplitude
        freq = 2.*pi.*f_reduced(B(1,1));                                                        % Find the frequency to the largest amplitude
        Omega(1,2) = freq;                                                                      % Set Omega2
    elseif(DYN.n_auto==2)
        [f,S,psi] = costarFFT(T,Z);                                                             % Do a fourier-transformation

        % Find first autonomous frequency
        ind_u1 = find(2.*pi.*f>0.95.*Omega(1,1),1);                                             % Find the frequency which is 5% lower than the estimated value for the autonomous frequency
        ind_o1 = find(2.*pi.*f>1.05.*Omega(1,1),1);                                             % Find the frequency which is 5% larger than the estimated value for the autonomous frequency
        f_reduced1 = f(ind_u1:ind_o1,1);                                                        % Make a window in the frequency domain
        S_reduced1 = S(ind_u1:ind_o1,1);
        [~,B1] = sort(S_reduced1,1,'descend');                                                  % Sort the frequencies in desceding order of amplitude
        freq1 = 2.*pi.*f_reduced1(B1(1,1));                                                     % Find the frequency to the largest amplitude
        Omega(1,1) = freq1;                                                                     % Set Omega1

        % Find second autonomous frequency
        ind_u2 = find(2.*pi.*f>0.95.*Omega(1,2),1);                                             % Find the frequency which is 5% lower than the estimated value for the autonomous frequency
        ind_o2 = find(2.*pi.*f>1.05.*Omega(1,2),1);                                             % Find the frequency which is 5% larger than the estimated value for the autonomous frequency
        f_reduced2 = f(ind_u2:ind_o2,1);                                                        % Make a window in the frequency domain
        S_reduced2 = S(ind_u2:ind_o2,1);
        [~,B2] = sort(S_reduced2,1,'descend');                                                  % Sort the frequencies in desceding order of amplitude
        freq2 = 2.*pi.*f_reduced2(B2(1,1));                                                     % Find the frequency to the largest amplitude
        Omega(1,2) = freq2;                                                                     % Set Omega2
    end

    %% Map data to grid in Hyper-Time
    %Find corrosponding torus coordinates
    theta1 = mod(Omega(1,1).*T,2*pi);                                                           % Define Hypertime coordinates
    theta2 = mod(Omega(1,2).*T,2*pi);                                                           % Define Hypertime coordinates

    %Define grid for interpolation
    X = linspace(0,2*pi,2*n_char);
    Y = linspace(0,2*pi,2*n_char);
    [XX,YY] = meshgrid(X,Y);
    
    %Interpolate data to grid -> TO DO: CHECK FOR DUPLICATE POINTS (WARNING BY MATLAB)
    for k=1:dim
        Cu = scatteredInterpolant(theta1,theta2,Z(:,k),'nearest');                              % Interpolate values from time-domain to hyper-time domain
        u(:,:,k) = Cu(XX,YY);
        clear Cu
    end

    %% Fit boundary values by a fourier-series
    %Fit Fourier series to boundary to find initial values for quasi-periodic shooting
    Xchar = linspace(0,2*pi*(1-1/(n_char+1)),n_char);
    for k=1:dim
        CU = fit(X.',u(:,1,k),'fourier8');                                                      % Do a fourier-fit of the boundary to generate intial values for the quasi-periodic shooting method
        U(:,k) = CU(Xchar).';
        clear CU
    end
    IV = reshape(U.',[dim*n_char,1]);                                                           % Reshape fitted data to get a vector

end

%% Save solution if autonomous frequencies are present
if(DYN.n_auto==0)
    obj.iv = IV;
elseif(DYN.n_auto==1)
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

elseif(DYN.n_auto==2)
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

%% Display information
disp('-------------------------------------------------------')
disp('---------------- Initial value found! -----------------')
disp('-------------------------------------------------------')

end


