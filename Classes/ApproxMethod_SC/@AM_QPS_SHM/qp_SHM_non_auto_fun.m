% Residual for Quasiperiodic Shooting - Heteronomous
%This function is a method of subclass AM_QPS_SHM. This function computes
%the residuum with QP-Shooting algorithm for first order nonautononomus ODE systems
%
%@obj:  ApproximationMethod subclass object
%@y:    solution curve vector
%@DYN:  DynamicalSystem class object
%
%@res:  residuum vector of evaluated ODE
%
function [f,J] = qp_SHM_non_auto_fun(obj,y,DYN)

%% Initialization
param = DYN.param;                                                                  % Get parameters
param{DYN.act_param} = y(end,1);                                                    % Set active parameter
Omega = DYN.non_auto_freq(y(end,1));                                                % Get frequencies
dim = DYN.dim;                                                                      % Get dimension of state-space
n_char = obj.n_char;                                                                % Get number of characteristics
x0 = reshape(y(1:end-1,1),[dim,n_char]);                                            % reshape solution vector

PHI(1,:) = obj.phi(1,:)./Omega(1,1);                                                % Rescale phi for multidimensional time
PHI(2,:) = obj.phi(2,:)./Omega(1,2);                                                % Rescale phi for multidimensional time

if(Omega(1,1) > Omega(1,2));  index = [1,2]; else;  index = [2,1]; end              % Define integration variable to minimize integration time
Ik = [0,2.*pi./Omega(1,index(1,1))];                                                % Set time span for integration

%% Time integration
dx = sqrt(eps)*(1 + max(abs(y(1:end-1,1)),[],1));                                   % Calculate differential to calculate Jacobian
INIT = repmat(reshape(x0,[dim,n_char]),[1,1,dim+1]);                                % Set initial values
for kk = 1:dim
    INIT(kk,:,kk+1) = x0(kk,:) + dx;                                                % Set initial values of perturbated characteristics
end

FW = @(t,z)obj.FcnWrapperODE5(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                    % Fcn Wrapper for time integration over all characteristics
[~,x] = obj.solver_function(FW,Ik,INIT,obj.odeOpts);                                % Time integration of all characteristics
xx = reshape(x(end,:).',[dim,n_char,dim+1]);                                        % Reshape end values of characteristics

G0 = repmat(reshape(xx(:,:,1),[dim*n_char,1]),[1,dim*n_char+1]);                    % Set values of G to end values of unperturbated characteristics
F = repmat(reshape(INIT(:,:,1),[dim*n_char,1]),[1,dim*n_char+1]);                   % Set F to unperturbated initial values
kjj = 0;
for lp = 1:dim:dim*n_char
    kjj = kjj+1;   
    G0(lp:lp+dim-1,lp+1:lp+dim) = permute(xx(:,kjj,2:dim+1),[1,3,2]);               % Set G to values of perturbated end values on mapped interval B(1,:)
    F(lp:lp+dim-1,lp+1:lp+dim) = permute(INIT(:,kjj,2:dim+1),[1,3,2]);              % Set F to values of perturbated initial values on initial interval [0,2pi(1-1/nchar)]
end

%% Remapping and interpolation
Theta = mod(obj.phi(index(1,2),:) + Ik(2)*Omega(1,index(1,2)),2*pi);                % Calculate end values of characteristics and map back to 0,2pi square
[B(1,:),I(1,:)] = sort(Theta);                                                      % Sort remapped values in ascending order

G01 = reshape(G0,[dim,n_char,dim*n_char+1]);                                        % Reshape G0 to be able to sort according to mapping in line 49
G = reshape(G01(:,I(1,:),:),[dim*n_char,dim*n_char+1]);                             % Map the end values of integration according to line 49

L0 = reshape(F.',[dim*(dim*n_char+1),n_char]);                                      % Reshape F to be able to interpolate in one step
C_3 = csape([obj.phi(index(1,2),:),2*pi],[L0,L0(:,1)],'periodic');                  % Interpolate initial values
H_temp = fnval(C_3,B(1,:));                                                         % Evaluate interpolated initial values on re-mapped points
H = reshape(H_temp,[n_char*dim+1,n_char*dim]).';                                    % Reshape back into original Form

%% Compute residuals and Jacobian
G1 = G-H;                                                                           % Set G1 to difference of start and end values (shooting residual)
f = G1(:,1);                                                                        % Set residual vector
J = (G1(:,2:end)-G1(:,1))./dx;                                                      % Calculate Jacobian

end



