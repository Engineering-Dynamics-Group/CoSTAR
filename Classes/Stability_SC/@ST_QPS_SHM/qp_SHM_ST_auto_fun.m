% Residual for Quasiperiodic Shooting - Full-autonomous
% This function is a method of subclass ST_QPS_SHM. 
% It computes the residuum with QP-Shooting algorithm for first order nonautononomus ODE systems
% Important: In contrast to the residuum function of AM_QPS_SHM the continuation parameter mu is not a variable
%
% @obj:  Stability subclass ST_QPS_SHM object
% @x:    Solution vector vector ( = y, but without mu)
% @mu:   Continuation parameter
% @DYN:  DynamicalSystem class object
%
% @f:    residuum vector of evaluated ODE
% @J:    Jacobian of f

function [f,J] = qp_SHM_ST_auto_fun(obj,x,mu,DYN)

%% Initialization
param = DYN.param;                                                                  % Get parameters
param{DYN.act_param} = mu;                                                          % Set active parameter
Omega(1,1) = x(end-1);                                                              % Get autonomous frequency 1 
Omega(1,2) = x(end);                                                                % Get autonomous frequency 2
dim = DYN.dim;                                                                      % Get dimension of state-space
n_char_st = obj.n_char_st;                                                          % Get number of characteristics
X = reshape(x(1:end-2),dim,n_char_st);                                              % reshape solution vector
AM_ST = obj.AM_ST;                                                                  % Object of AM_QPS_SHM

f = zeros(dim*n_char_st+2,1);                                                       % Initialize residual vector
J = zeros(dim*n_char_st+2,dim*n_char_st+2);                                         % Initialize Jacobian

f1 = AM_ST.Y_old{1,3};                                                              % Gradient of reference solution by theta2
h1 = AM_ST.Y_old{1,4};                                                              % Gradient of reference solution by theta1 
reso_phase = AM_ST.reso_phase;                                                      % Time resolution of phase condition

phi = AM_ST.phi;                                                                    % Get phi
PHI(1,:) = phi(1,:)./Omega(1,1);                                                    % Rescale phi for multidimensional time
PHI(2,:) = phi(2,:)./Omega(1,2);                                                    % Rescale phi for multidimensional time

Ik = linspace(0,2.*pi./Omega(1,1),reso_phase);                                      % Set time span for integration

%% Time integration
dx = sqrt(eps)*(1 + max(abs(x),[],1));                                              % Define differential for Jacobian
INIT = repmat(reshape(X,[dim,n_char_st]),[1,1,dim+1]);                              % Set initial values
for kk = 1:dim
    INIT(kk,:,kk+1) = X(kk,:) + dx;                                                 % Set initial values of perturbated characteristics
end

FW = @(t,z)AM_ST.FcnWrapperODE5(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                  % Fcn Wrapper for time integration over all characteristics
[~,z] = AM_ST.solver_function(FW,Ik,INIT,AM_ST.odeOpts);                            % Time integration of all characteristics
zz = reshape(z(end,:).',[dim,n_char_st,dim+1]);                                     % Reshape end values of characteristics
z_sort = reshape(z,[reso_phase,dim,n_char_st,dim+1]);                               % Save data of time integration for phase condition

G0 = repmat(reshape(zz(:,:,1),[dim*n_char_st,1]),[1,dim*n_char_st+1]);              % Set values of G to end values of unperturbated characteristics
F = repmat(reshape(INIT(:,:,1),[dim*n_char_st,1]),[1,dim*n_char_st+1]);             % Set F to unperturbated initial values
kjj = 0;
for lp = 1:dim:dim*n_char_st
    kjj = kjj+1;   
    G0(lp:lp+dim-1,lp+1:lp+dim) = permute(zz(:,kjj,2:dim+1),[1,3,2]);               % Set G to values of perturbated end values on mapped interval B(1,:)
    F(lp:lp+dim-1,lp+1:lp+dim) = permute(INIT(:,kjj,2:dim+1),[1,3,2]);              % Set F to values of perturbated initial values on initial interval [0,2pi(1-1/nchar)]
end

%% Remapping and interpolation
Theta = mod(phi(2,:) + 2*pi*Omega(1,2)/Omega(1,1),2*pi);                            % Calculate end values of characteristics and map back to 0,2pi square
[B(1,:),I(1,:)] = sort(Theta);                                                      % Sort remapped characteristics in ascending order

G01 = reshape(G0,[dim,n_char_st,dim*n_char_st+1]);                                  % Reshape G0 to be able to sort according to mapping in line 57
G = reshape(G01(:,I(1,:),:),[dim*n_char_st,dim*n_char_st+1]);                       % Map the end values of integration according to line 57
    
L0 = reshape(F.',[dim*(dim*n_char_st+1),n_char_st]);                                % Reshape F to be able to interpolate in one step
C_3 = csape([phi(2,:),2*pi],[L0,L0(:,1)],'periodic');                               % Interpolate initial values
H_temp = fnval(C_3,B(1,:));                                                         % Evaluate interpolated initial values on re-mapped points
H = reshape(H_temp,[n_char_st*dim+1,n_char_st*dim]).';                              % Reshape back into original Form

f2 = permute(z_sort(:,:,:,1),[1,3,2,4]);                                            % Fetch solution with unperturbated initial values for phase condition
f(end-1,1) = AM_ST.poincare_int(f2,h1,Omega,Ik);                                      % Set phase-condition for autonomous frequency 1
f(end,1) = AM_ST.poincare_int(f2,f1,Omega,Ik);                                        % Set phase-condition for autonomous frequency 2

k_1 = 1;
for jk = 1:n_char_st
    for jjk = 1:dim
        f3 = f2;                                                                    % Set values to unperturbated values
        f3(:,jk,:) = z_sort(:,:,jk,1+jjk);                                          % Perturbate value of one characteristic
        P1 = AM_ST.poincare_int(f3,h1,Omega,Ik);                                      % Evaluate integral poincare phase condition for perturbated initial conditions of autonoumous frequency 1
        J(end-1,k_1) = (P1-f(end-1,1))/dx;                                          % Caluclate Jacobian dP1/dx 
        
        P2 = AM_ST.poincare_int(f3,f1,Omega,Ik);                                      % Evaluate integral poincare phase condition for perturbated initial conditions of autonoumous frequency 2
        J(end,k_1) = (P2-f(end,1))/dx;                                              % Caluclate Jacobian dP2/dx 
        k_1 = k_1 + 1;
    end
end

%% Compute residuals and Jacobian
G1 = G-H;                                                                           % Set G1 to difference of start and end values (shooting residual)
f(1:end-2,1) = G1(:,1);                                                             % Set residual vector
J(1:end-2,1:end-2) = (G1(:,2:end)-G1(:,1))./dx;                                     % Calculate Jacobian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Autonomous frequency 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
Om1_old = x(end-1);                                                                 % Reset first autonomous frequency to unpertubated value
x(end-1) = x(end-1)+dx;                                                             % Perurbate first autonomous frequency
Omega = [x(end-1),x(end)];                                                          % Get frequencies

Theta1 = mod(phi(2,:) + 2*pi*Omega(1,2)/Omega(1,1),2*pi);                           % Calculate end values of characteristics and map back to 0,2pi square for perturbated omega1
[B1,I1] = sort(Theta1);                                                             % Sort remapped characteristics in ascending order

Ik1 = linspace(0,2.*pi./Omega(1,1),reso_phase);                                     % Set integration interval for perturbated Omega1
PHI(1,:) = phi(1,:)./Omega(1,1);                                                    % Set phase difference of characteristics
PHI(2,:) = phi(2,:)./Omega(1,2);                                                

%% Time integration
FW = @(t,z)AM_ST.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                  % Fcn Wrapper for time integration over all characteristics
[~,z1] = AM_ST.solver_function(FW,Ik1,INIT(:,:,1),AM_ST.odeOpts);                   % Time integration of all characteristics
z1_sort = reshape(z1,[reso_phase,dim,n_char_st]);                                   % Save data of time integration for phase condition
G0temp = reshape(z1(end,:).',[dim,n_char_st]);                                      % Reshape end values of time integration to be able to sort them                                     
G02 = reshape(G0temp(:,I1),[dim*n_char_st,1]);                                      % Sort the end values according to mapping

f4 = permute(z1_sort,[1,3,2]);                                                      % Perturbed to [reso_phase,n_char,dim,perturb]
L1 = reshape(F(:,1).',[dim,n_char_st]);                                             % Reshape F to be able to interpolate in one step
C_01 = csape([phi(2,:),2*pi],[L1,L1(:,1)],'periodic');                              % Interpolate initial values
H01 = reshape(fnval(C_01,B1(1,:)),[dim*n_char_st,1]);                               % Evaluate interpolated initial values on re-mapped points

J(end-1,end-1) = (AM_ST.poincare_int(f4,h1,Omega,Ik1)-f(end-1,1))./dx;              % Set values of of Jacobian dP1/domega1
J(end,end-1) = (AM_ST.poincare_int(f4,f1,Omega,Ik1)-f(end,1))./dx;                  % Set values of last of Jacobian dP2/domega1             

G2 = G02 - H01(:,1);                                                                % Set residual for perturbate autonomous frequency 1
J(1:dim*n_char_st,end-1) = (G2-G1(:,1))./dx;                                        % Set values of second to last column of Jacobian df/domega1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Autonomous frequency 2%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(end-1) = Om1_old;                                                                 % Reset first autonomous frequency
x(end) = x(end)+dx;                                                                 % Perurbate second autonomous frequency
Omega = [x(end-1),x(end)];                                                          % Set frequencies to perturbated frequency

Theta2 = mod(phi(2,:) + 2*pi*Omega(1,2)/Omega(1,1),2*pi);                           % Calculate end values of characteristics and map back to 0,2pi square for perturbated omega2
[B2,I2] = sort(Theta2);                                                             % Sort mapped values for characteristics

Ik2 = linspace(0,2.*pi./Omega(1,1),reso_phase);                                     % Define integation interval
PHI(1,:) = phi(1,:)./Omega(1,1);                                                    % Set phase difference of characteristics
PHI(2,:) = phi(2,:)./Omega(1,2);

%% Time integration
FW = @(t,z)AM_ST.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                  % Fcn Wrapper for time integration over all characteristics
[~,z2] = AM_ST.solver_function(FW,Ik2,INIT(:,:,1),AM_ST.odeOpts);                   % Time integration of all characteristics
z2_sort = reshape(z2,[reso_phase,dim,n_char_st]);                                   % Save data of time integration for phase condition

f5 = permute(z2_sort,[1,3,2]);                                                      % Perturbed to [reso_phase,n_char,dim,perturb]
L2 = reshape(F(:,1).',[dim,n_char_st]);                                             % Reshape F to be able to interpolate in one step
C_02 = csape([phi(2,:),2*pi],[L2,L2(:,1)],'periodic');                              % Interpolate initial values
H02 = reshape(fnval(C_02,B2(1,:)),[dim*n_char_st,1]);                               % Evaluate interpolated initial values on re-mapped points

J(end-1,end) = (AM_ST.poincare_int(f5,h1,Omega,Ik2)-f(end-1,1))./dx;                % Set values of of Jacobian dP1/domega2
J(end,end) = (AM_ST.poincare_int(f5,f1,Omega,Ik2)-f(end,1))./dx;                    % Set values of of Jacobian dP2/domega2

G03temp = reshape(z2(end,:).',[dim,n_char_st]);                                     % Reshape end values of time integration to be able to sort them          
G03 = reshape(G03temp(:,I2),[dim*n_char_st,1]);                                     % Sort the end values according to mapping

G3 = G03 - H02(:,1);                                                                % Set residual for perturbated autonomous frequency 2
J(1:dim*n_char_st,end) = (G3-G1(:,1))./dx;                                          % Set values of last column of Jacobian df/domega2

end