% Reference solution for the mixed multiple shooting phase condition
% This function is a method of subclass AM_QPS_SHM.
% Integrates each reference segment from its own shooting node and stores
% the derivative with respect to the free phase theta2 on the same grid.
%
%@obj:     ApproximationMethod subclass object
%@Z0_nodes:Shooting nodes of size [dim,n_shoot,n_char]
%@Omega:   Reference frequency vector [Omega1,Omega2]
%@param:   Parameters of the reference solution
%@DYN:     DynamicalSystem class object
%
function obj = updatePhaseReference(obj,Z0_nodes,Omega,param,DYN)

dim = DYN.dim;                                                              % Get dimension of state-space
n_char = obj.n_char;                                                        % Get number of characteristics
T = 2*pi/Omega(1,1);                                                        % Integration period for the prescribed frequency
[Ik_phase,tau_phase] = obj.getPhaseGrids(T);                                % Generate the same grids as the mixed residual
PHI = obj.phi./Omega(:);                                                    % Rescale phi for multidimensional time
FW = @(t,z)obj.FcnWrapperODE2(t,z,@(t,z)DYN.rhs(t,z,param),PHI);
W_phase = cell(1,obj.n_shoot);
F1_phase = cell(1,obj.n_shoot);
tau = linspace(0,1,obj.reso_phase);
W = zeros(obj.reso_phase,n_char,dim);
F1 = zeros(obj.reso_phase,n_char,dim);

for i = 1:obj.n_shoot
    IV = reshape(Z0_nodes(:,i,:),[dim*n_char,1]);                           % Fetch initial conditions of the current reference segment
    [~,V] = obj.solver_function(FW,Ik_phase{i},IV,obj.odeOpts);             % Time integration of all characteristics over the current interval
    W_i = permute(reshape(V,[numel(Ik_phase{i}),dim,n_char]),[1,3,2]);      % Reshape solution
    F1_i = zeros(size(W_i));                                                % Preallocate the matrix that will be filled with the derivatives of the reference solution W_i
    for j = 1:dim
        F1_i(:,:,j) = gradient(W_i(:,:,j),obj.phi(2,:),tau_phase{i});       % Calculate derivative of reference solution with respect to \theta_2 
    end
    
    W_phase{i} = W_i;                                                       % Save reference solution on the current interval, i.e.
    F1_phase{i} = F1_i;                                                     % Save derivative on the same normalized time grid

    % Retain the numeric Y_old representation.
    rows = find(tau>=tau_phase{i}(1) & tau<=tau_phase{i}(end));
    W(rows,:,:) = reshape(interp1(tau_phase{i}, ...
        reshape(W_i,size(W_i,1),[]),tau(rows)),[numel(rows),n_char,dim]);
    F1(rows,:,:) = reshape(interp1(tau_phase{i}, ...
        reshape(F1_i,size(F1_i,1),[]),tau(rows)),[numel(rows),n_char,dim]);
end

ref.tau = tau_phase;
ref.values = W_phase;
ref.gradient = F1_phase;
obj.phase_reference = ref;                                               % Replace reference data only between nonlinear solution processes
obj.Y_old{1,1} = obj.phi(2,:);                                            % Save integration interval
obj.Y_old{1,2} = W;                                                       % Save reference solution
obj.Y_old{1,3} = F1;                                                      % Save derivative of reference solution with respect to theta2
end
