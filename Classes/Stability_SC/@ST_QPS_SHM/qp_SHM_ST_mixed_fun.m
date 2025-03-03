% Residual for Quasiperiodic Shooting - Mixed Case
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

function [f,J] = qp_SHM_ST_mixed_fun(obj,x,mu,DYN)

%% Initialization
param = DYN.param;                                                                  % Get parameters
param{DYN.act_param} = mu;                                                          % Set active parameter
Omega(1,1) = DYN.non_auto_freq(mu);                                                 % Get non-autonomous frequency
Omega(1,2) = x(end);                                                                % Get autonomous frequency
dim = DYN.dim;                                                                      % Get dimension of state-space
n_char_st = obj.n_char_st;                                                          % Get number of characteristics
X = reshape(x(1:end-1),dim,n_char_st);                                              % reshape solution vector
AM_ST = obj.AM_ST;                                                                  % Object of AM_QPS_SHM

f = zeros(dim*n_char_st+1,1);                                                       % Initialize residual vector
J = zeros(dim*n_char_st+1);                                                         % Initialize Jacobian

F1 = AM_ST.Y_old{1,3};                                                              % Gradient of reference solution
reso_phase = AM_ST.reso_phase;                                                      % Time resolution of phase condition

phi = AM_ST.phi;                                                                    % Get phi
PHI(1,:) = phi(1,:)./Omega(1,1);                                                    % Rescale phi for multidimensional time
PHI(2,:) = phi(2,:)./Omega(1,2);                                                    % Rescale phi for multidimensional time

Ik = linspace(0,2.*pi./Omega(1,1),reso_phase);                                      % Set time span for integration

%% Time integration
dx = sqrt(eps)*(1 + max(abs(x),[],1));
INIT = repmat(reshape(X,[dim,n_char_st]),[1,1,dim+1]);                              % Set initial values
for kk = 1:dim
    INIT(kk,:,kk+1) = X(kk,:) + dx;                                                 % Set initial values of perturbated characteristics
end

FW = @(t,z)AM_ST.FcnWrapperODE5(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                  % Fcn Wrapper for time integration over all characteristics
[~,x] = AM_ST.solver_function(FW,Ik,INIT,AM_ST.odeOpts);                            % Time integration of all characteristics
xx = reshape(x(end,:).',[dim,n_char_st,dim+1]);                                     % Reshape end values of characteristics
x_sort = reshape(x,[reso_phase,dim,n_char_st,dim+1]);

G0 = repmat(reshape(xx(:,:,1),[dim*n_char_st,1]),[1,dim*n_char_st+1]);              % Set values of G to end values of unperturbated characteristics
F = repmat(reshape(INIT(:,:,1),[dim*n_char_st,1]),[1,dim*n_char_st+1]);             % Set F to unperturbated initial values
kjj = 0;
for lp = 1:dim:dim*n_char_st
    kjj = kjj+1;   
    G0(lp:lp+dim-1,lp+1:lp+dim) = permute(xx(:,kjj,2:dim+1),[1,3,2]);               % Set G to values of perturbated end values on mapped interval B(1,:)
    F(lp:lp+dim-1,lp+1:lp+dim) = permute(INIT(:,kjj,2:dim+1),[1,3,2]);              % Set F to values of perturbated initial values on initial interval [0,2pi(1-1/nchar)]
end

%% Remapping and interpolation
Theta = mod(phi(2,:) + 2*pi*Omega(1,2)/Omega(1,1),2*pi);                            % Calculate end values of characteristics and map back to 0,2pi square
[B(1,:),I(1,:)] = sort(Theta);                                                      % Sort remapped values in ascending order

G01 = reshape(G0,[dim,n_char_st,dim*n_char_st+1]);                                  % Reshape G0 to be able to sort according to mapping in line 57
G = reshape(G01(:,I(1,:),:),[dim*n_char_st,dim*n_char_st+1]);                       % Map the end values of integration according to line 57

L0 = reshape(F.',[dim*(dim*n_char_st+1),n_char_st]);                                % Reshape F to be able to interpolate in one step
C_3 = csape([phi(2,:),2*pi],[L0,L0(:,1)],'periodic');                               % Interpolate initial values
H_temp = fnval(C_3,B(1,:));                                                         % Evaluate interpolated initial values on re-mapped points
H = reshape(H_temp,[n_char_st*dim+1,n_char_st*dim]).';                              % Reshape back into original Form

%% Phase
% Time delta & distorted right boundry                                                                            
Omega_22 = Omega(1,2) + dx;                                                         % Perturbated autonomous omega

Theta = mod(phi(2,:)+2*pi*Omega_22/Omega(1,1),2*pi);                                % Mapping of endvalues of theta for perturbated autonomous frequency
[B1,~] = sort(Theta);                                                               % Sorting according to mapping

L1 = reshape(F(:,1).',[dim,n_char_st]);                                             % Reshape F to be able to interpolate in one step
C_4 = csape([phi(2,:),2*pi],[L1,L1(:,1)],'periodic');                               % Interpolate initial values
H_temp1 = fnval(C_4,B1(1,:));                                                       % Evaluate interpolated initial values on re-mapped points
H1 = reshape(H_temp1,[n_char_st*dim,1]);                                            % Reshape back into original Form

F2 = permute(x_sort(:,:,:,1),[1,3,2,4]);                                            % Fetch solution with unperturbated initial values for phase condition
f(end,1) = AM_ST.poincare_int(F2,F1,Omega,Ik);                                      % Set phase-condition for autonomous frequency

k_1 = 1;
for jk = 1:n_char_st
    for jjk = 1:dim
        F3 = F2;                                                                    % Set values to unperturbated values
        F3(:,jk,:) = x_sort(:,:,jk,1+jjk);                                          % Perturbate value of one characteristic
        P1 = AM_ST.poincare_int(F3,F1,Omega,Ik);                                    % Evaluate integral poincare phase condition for perturbated initial conditions
        J(end,k_1) = (P1-f(end,1))/dx;                                              % Caluclate Jacobian dP/dx 
        k_1 = k_1 + 1;
    end
end

%% Compute residuals and Jacobian
G1 = [G,G(:,1)]-[H,H1];                                                             % Set G1 to difference of start and end values (shooting residual)
f(1:end-1,1) = G1(:,1);                                                             % Set residual vector
J(1:end-1,:) = sparse((G1(:,2:end)-G1(:,1))./dx);                                   % Calculate Jacobian

end