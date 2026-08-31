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
n_shoot = obj.n_shoot;                                                              % Get number of shooting points
n_state = dim*n_shoot*n_char;

%% Switch old/new single shooting residual
old_single_shooting = false;


Z0 = reshape(y(1:n_state),[dim,n_shoot,n_char]);                         % reshape solution vector and exclude the active parameter [located at y(n_state+1)]


PHI(1,:) = obj.phi(1,:)./Omega(1,1);                                                % Rescale phi for multidimensional time
PHI(2,:) = obj.phi(2,:)./Omega(1,2);                                                % Rescale phi for multidimensional time

if(Omega(1,1) > Omega(1,2));  index = [1,2]; else;  index = [2,1]; end              % Define integration variable to minimize integration time
Ik = [0,2.*pi./Omega(1,index(1,1))];                                                % Set time span for integration


%% single shooting with old calculation of the Jacobian
if old_single_shooting && n_shoot == 1
    %% Time integration
    dx = sqrt(eps)*(1 + max(abs(y(1:end-1,1)),[],1));                                   % Calculate differential to calculate Jacobian
    INIT = repmat(reshape(Z0,[dim,n_char]),[1,1,dim+1]);                                % Set initial values
    for kk = 1:dim
        INIT(kk,:,kk+1) = Z0(kk,:) + dx;                                                % Set initial values of perturbated characteristics
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
        F(lp:lp+dim-1,lp+1:lp+dim) = permute(INIT(:,kjj,2:dim+1),[1,3,2]);              % Set F to values of perturbated initial values on initial interval [0,2pi(1-1/n_char)]
    end

    %% Remapping and interpolation
    Theta = mod(obj.phi(index(1,2),:) + Ik(2)*Omega(1,index(1,2)),2*pi);                % Calculate end values of characteristics and map back to 0,2pi square
    [sortedPhases(1,:),sortindex(1,:)] = sort(Theta);                                                      % Sort remapped values in ascending order
    
    G01 = reshape(G0,[dim,n_char,dim*n_char+1]);                                        % Reshape G0 to be able to sort according to mapping in line 49
    G = reshape(G01(:,sortindex(1,:),:),[dim*n_char,dim*n_char+1]);                             % Map the end values of integration according to line 49
    
    L0 = reshape(F.',[dim*(dim*n_char+1),n_char]);                                      % Reshape F to be able to interpolate in one step
    C_3 = csape([obj.phi(index(1,2),:),2*pi],[L0,L0(:,1)],'periodic');                  % Interpolate initial values
    H_temp = fnval(C_3,sortedPhases(1,:));                                                         % Evaluate interpolated initial values on re-mapped points
    H = reshape(H_temp,[n_char*dim+1,n_char*dim]).';                                    % Reshape back into original Form
    
    %% Compute residuals and Jacobian
    G1 = G-H;                                                                           % Set G1 to difference of start and end values (shooting residual)
    f = G1(:,1);                                                                        % Set residual vector
    J = sparse((G1(:,2:end)-G1(:,1))./dx);                                                      % Calculate Jacobian



    
%% multiple shooting / new single shooting
else
    % Calculate residual and Jacobian for the unperturbed shooting-node vector
    [f,J] = computeMSResidual(obj,y,DYN);

end

end

function [R,J] = computeMSResidual(obj,y,DYN)
%COMPUTEMSRESIDUAL
% Computes only the multiple-shooting residual for a non-autonomous QPS.

%% Initialization

% Extract needed parameters from the DYN object
param = DYN.param;
param{DYN.act_param} = y(end,1);
Omega = DYN.non_auto_freq(y(end,1));
dim = DYN.dim;
n_char = obj.n_char;
n_shoot = obj.n_shoot;
n_state = dim*n_shoot*n_char;

Z0_nodes = reshape(y(1:n_state),[dim,n_shoot,n_char]);                 % Reshape the initial value vector into the corresponding matrix

PHI = zeros(2,n_char);
PHI(1,:) = obj.phi(1,:)./Omega(1,1);
PHI(2,:) = obj.phi(2,:)./Omega(1,2);

if Omega(1,1) > Omega(1,2)
    index = [1,2];
else
    index = [2,1];
end
T = 2*pi/Omega(1,index(1));
Ik = [0,T];


%% perturbe initial values
dz = sqrt(eps)*(1 + max(max(max(abs(Z0_nodes)))));                                % Define a a small perturbation value

Z0_perturbed = repmat(Z0_nodes, [1,1,1, dim+1]);                                       % Expand Z0 to next perturbe each state seperately

% For last index = 1 we retain the unperturbed initial value. 
% For last index > 1 we get a pertubation in the state index - 1.
for i = 1:dim
    Z0_perturbed(i, :, :, i+1) = Z0_nodes(i, :, :) + dz;                              
end

%% Sub-intervals

dT = (Ik(2)-Ik(1))/n_shoot;                                                 % Define the length of one multiple shooting interval

T_start = Ik(1) + (0:n_shoot-1)*dT;                                         % Define initial start value for forward integration
T_end   = Ik(1) + (1:n_shoot)*dT;                                           % Define initial end value for forward integration

Z_SubEnd = zeros(dim,n_shoot,n_char, dim+1,'like',Z0_perturbed);            % Initialising Z_SubEnd with 'like' ensures the matrix has the same data type, sparsity and attributes of Z0_perturbed

FW = @(t,z)obj.FcnWrapperODE5(t,z,@(t,z)DYN.rhs(t,z,param),PHI);

%% Integrate all shooting intervals



for i = 1:n_shoot                             

    Z_start = reshape(Z0_perturbed(:,i,:,:), ...                            % Reshape the initial values used for this shooting interval into a  vector
        [dim*n_char*(dim+1),1]);                                            % as this is more natural for ODEs solvers.
        

    [~,Z] = obj.solver_function( ...
        FW, ...
        [T_start(i),T_end(i)], ...
        Z_start, ...
        obj.odeOpts);

    Z_SubEnd(:,i,:,:) = reshape(Z(end,:),[dim,1,n_char,dim+1]);             % Reshape the solution given as a vector into a matrix

end

% Z_end(:, :, :) = Z_SubEnd(:, end, :, :);
Z_end = reshape(Z_SubEnd(:,end,:,:),[dim,n_char,dim+1]);


% Get final values of each characteristic 
finPhases = mod(obj.phi(index(2), :) + T*Omega(index(2)), 2*pi);

% Sort the end phases and get permutation matrix to reorder the states
% accordingly
[finPhasesSorted, sortindex] = sort(finPhases);
inverseSortindex(sortindex) = 1:numel(sortindex);                     % Get the inverse of the sorting index for later use

% Arrange integrated end values for the boundary residual
if n_shoot == 1
    Z_endDiagMatr = obj.perturbMatr(Z_end,dim, n_char);
else
    TempMatr = reshape(Z_end, [], dim + 1);
    Z_endDiagMatr = repmat(TempMatr(:, 1), [1, dim*n_char+1]);
end

Z0DiagMatr = obj.perturbMatr(permute(Z0_perturbed(:, 1, :, :), [1, 3, 4, 2]), dim, n_char);
Z0Vec = reshape(Z0DiagMatr.', [dim*(dim*n_char + 1), n_char]);

% Reshape the diagonally perturbed matrix
Z_endReshaped = reshape(Z_endDiagMatr, [dim, n_char, dim*n_char+1]);
Z_endSorted = reshape(Z_endReshaped(:, sortindex, :), [dim*n_char, dim*n_char+1]);


%% periodic boundary conditions
x_values = [obj.phi(index(2), :), 2*pi];                                  % Append 2*pi
y_values = [Z0Vec, Z0Vec(:, 1)];                                            % Append first value

spline = csape(x_values, y_values, 'periodic');                             

evalSpline = fnval(spline, finPhasesSorted);                                % Evaluate the interpolated polynomial of the initial states at the sorted final phases
evalSplineMatr = reshape(evalSpline, [n_char*dim + 1, n_char*dim]).';         % Reshape the vector into a matrix

ResidualMatr = Z_endSorted - evalSplineMatr;                          % Define residual of the end points of the characteristics

ResidualMatrReshaped = reshape(ResidualMatr, [dim, n_char, dim*n_char+1]);  % Reshape the Residual  
ResidualMatr = reshape( ResidualMatrReshaped(:, inverseSortindex, :), ... % Reverse the sorting
    [dim*n_char,dim*n_char+1]);


%% Assemble the complete residual matrix

% Indices corresponding to Z0_nodes(dim,n_shoot,n_char)
node_index = reshape(1:n_state,[dim,n_shoot,n_char]);

% Rows containing the quasiperiodic boundary equations
boundary_rows = reshape(node_index(:,end,:),[],1);

% Embed the boundary residual into the full multiple-shooting matrix
expandedResidualMatr = obj.expandResidualMatr( ...
    ResidualMatr,dim,n_char,n_shoot);

if n_shoot == 1

    % For single shooting, ResidualMatr already contains the flow
    % derivatives and the spline derivatives.
    FullResidualMatr = expandedResidualMatr;

else

    % Rows containing the inner continuity equations
    interior_rows = reshape( ...
        node_index(:,1:end-1,:),[],1);

    % Corresponding rows of the following shooting nodes
    next_node_rows = reshape( ...
        node_index(:,2:end,:),[],1);

    % Matrices containing the perturbations of all shooting nodes
    Z0SubDiagMatr = obj.perturbMatrSub( ...
        Z0_perturbed,dim,n_char,n_shoot);

    Z_SubEndDiagMatr = obj.perturbMatrSub( ...
        Z_SubEnd,dim,n_char,n_shoot);

    FullResidualMatr = zeros( ...
        n_state,n_state+1,'like',ResidualMatr);

    %% Inner continuity equations
    % Z_end(i) - Z_start(i+1) = 0
    FullResidualMatr(interior_rows,:) = ...
        Z_SubEndDiagMatr(interior_rows,:) ...
        - Z0SubDiagMatr(next_node_rows,:);

    %% Quasiperiodic boundary equations
    % This part contains the spline derivatives with respect to the first
    % shooting nodes.
    FullResidualMatr(boundary_rows,:) = ...
        expandedResidualMatr(boundary_rows,:);

    % Add only the variations of the integrated endpoints with respect
    % to the last shooting nodes. Subtracting the first column removes
    % the unperturbed endpoint values.
    boundaryFlowDelta = ...
        Z_SubEndDiagMatr(boundary_rows,:) ...
        - repmat( ...
            Z_SubEndDiagMatr(boundary_rows,1), ...
            [1,n_state+1]);

    FullResidualMatr(boundary_rows,:) = ...
        FullResidualMatr(boundary_rows,:) ...
        + boundaryFlowDelta;

end

%% Assemble residual and Jacobian

R = FullResidualMatr(:,1);
R = R(:);

J = sparse( ...
    (FullResidualMatr(:,2:end) ...
    - FullResidualMatr(:,1))./dz);

end
