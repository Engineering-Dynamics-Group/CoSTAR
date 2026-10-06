% Residual for Quasiperiodic Shooting - Mixed Case
%This function is a method of subclass AM_QPS_SHM. This function computes
%the residuum with QP-Shooting algorithm for first order mixed ODE systems
%with one non-autonomous and one autonomous frequency
%
%@obj:  ApproximationMethod subclass object
%@y:    solution curve vector containing auto frequency
%@DYN:  DynamicalSystem class object
%
%@res:  residuum vector of evaluated ODE
%



function [f,J] = qp_SHM_mixed_fun(obj,y,DYN)
%% Initialization
% Extract needed parameters from the DYN object
param = DYN.param;                                                                % Get parameters
param{DYN.act_param} = y(end,1);                                                  % Set active parameter
Omega = [DYN.non_auto_freq(y(end,1)),y(end-1,1)];
dim = DYN.dim;                                                                    % Get dimension of state-space
n_char = obj.n_char;                                                              % Get number of characteristics
n_shoot = obj.n_shoot;                                                            % Get number of shooting points
n_state = dim*n_shoot*n_char;

Z0_nodes = reshape(y(1:n_state),[dim,n_shoot,n_char]);                      % Reshape the initial value vector into the corresponding matrix

PHI = zeros(2,n_char);
PHI(1,:) = obj.phi(1,:)./Omega(1,1);                                              % Rescale phi_1 w.r.t. \theta_1
PHI(2,:) = obj.phi(2,:)./Omega(1,2);                                              % Rescale phi_2 w.r.t. \theta_2

% Integrate over the prescribed period; the free phase is \theta_2.
T = 2*pi/Omega(1,1);
[Ik_phase,~,F1_phase] = obj.getPhaseGrids(T);                                   % Set time grids and interpolate the fixed reference gradient
F_phase = cell(1,n_shoot);                                                % Store solution values for the phase condition

%% Residual-only evaluation
if nargout < 2


    % Only unperturbed trajectories are needed.
    Z_SubEnd = zeros(dim,n_shoot,n_char,'like',Z0_nodes);

    FW = @(t,z)obj.FcnWrapperODE2( ...
        t,z,@(t,z)DYN.rhs(t,z,param),PHI);

    for i = 1:n_shoot

        Z_start = reshape( ...
            Z0_nodes(:,i,:),[dim*n_char,1]);

        [~,Z] = obj.solver_function( ...
            FW, ...
            Ik_phase{i}, ...
            Z_start, ...
            obj.odeOpts);

        Z_SubEnd(:,i,:) = reshape( ...
            Z(end,:).',[dim,1,n_char]);
        F_phase{i} = permute(reshape(Z,[numel(Ik_phase{i}),dim,n_char]), ...
            [1,3,2]);                                                    % Save values along the current characteristic segments

    end

    %% Quasiperiodic boundary residual
    Z_end = reshape( ...
        Z_SubEnd(:,end,:),[dim,n_char]);

    finPhases = mod( ...
        obj.phi(2,:) + T*Omega(2), ...
        2*pi);

    [finPhasesSorted,sortindex] = sort(finPhases);

    inverseSortindex = zeros(size(sortindex));
    inverseSortindex(sortindex) = 1:n_char;

    Z_boundary = reshape( ...
        Z0_nodes(:,1,:),[dim,n_char]);

    boundarySpline = csape( ...
        [obj.phi(2,:),2*pi], ...
        [Z_boundary,Z_boundary(:,1)], ...
        'periodic');

    Z_boundaryMapped = reshape( ...
        fnval(boundarySpline,finPhasesSorted), ...
        [dim,n_char]);

    boundaryResidualSorted = ...
        Z_end(:,sortindex) - Z_boundaryMapped;

    boundaryResidual = ...
        boundaryResidualSorted(:,inverseSortindex);

    %% Assemble complete residual vector
    node_index = reshape( ...
        1:n_state,[dim,n_shoot,n_char]);

    boundary_rows = reshape( ...
        node_index(:,end,:),[],1);

    f = zeros(n_state,1,'like',Z0_nodes);

    if n_shoot > 1

        interior_start_nodes = reshape( ...
            node_index(:,1:end-1,:),[],1);

        interior_end_nodes = reshape( ...
            node_index(:,2:end,:),[],1);

        f(interior_start_nodes) = ...
            Z_SubEnd(interior_start_nodes) ...
            - Z0_nodes(interior_end_nodes);

    end

    f(boundary_rows) = boundaryResidual(:);

    f(end+1,1) = obj.poincare_int_multiple(F_phase,F1_phase,Omega,Ik_phase);

    return
end

%% Compute Residual and Jacobian

%% Perturb initial values
dz = sqrt(eps)*(1 + max(max(max(abs(Z0_nodes)))));                          % Define a small perturbation value

Z0_perturbed = repmat(Z0_nodes, [1,1,1, dim+1]);                            % Z0_perturbed contains the unperturbed initial values and a new dimension for later perturbations.

% For last index = 1 we retain the unperturbed initial value.
% For last index > 1 we get a pertubation in the state index - 1.
for i = 1:dim
    Z0_perturbed(i, :, :, i+1) = Z0_nodes(i, :, :) + dz;
end

%% Sub-intervals
F_perturbed = cell(1,n_shoot);                                            % Store perturbed trajectories for the phase derivatives

Z_SubEnd = zeros(dim,n_shoot,n_char, dim+1,'like',Z0_perturbed);            % Initializing Z_SubEnd with 'like' ensures the matrix has the same data type, sparsity and attributes of Z0_perturbed

FW = @(t,z)obj.FcnWrapperODE5(t,z,@(t,z)DYN.rhs(t,z,param),PHI);                  % Fcn Wrapper for time integration over all characteristics

%% Integrate all shooting intervals
for i = 1:n_shoot
    Z_start = reshape(Z0_perturbed(:,i,:,:), ...                            % Reshape the initial values used for this shooting interval into a  vector
        [dim*n_char*(dim+1),1]);                                            % since this representation is more natural for ODEs solvers.


    [~,Z] = obj.solver_function( ...
        FW, ...
        Ik_phase{i}, ...
        Z_start, ...
        obj.odeOpts);

    F_perturbed{i} = reshape(Z,[numel(Ik_phase{i}),dim,n_char,dim+1]);
    F_phase{i} = permute(F_perturbed{i}(:,:,:,1),[1,3,2]);                  % Fetch solution with unperturbed initial values for phase condition
    Z_SubEnd(:,i,:,:) = reshape(Z(end,:),[dim,1,n_char,dim+1]);             % Reshape the solution given as a vector into a matrix
end

% Z_end(:, :, :) = Z_SubEnd(:, end, :, :);
Z_end = reshape(Z_SubEnd(:,end,:,:),[dim,n_char,dim+1]);

finPhases = mod(obj.phi(2, :) + T*Omega(2), 2*pi);            % Calculate the mapped final phase of each characteristic.

[finPhasesSorted, sortindex] = sort(finPhases);                             % Sort the end phases and get index vector to reorder the states accordingly
inverseSortindex(sortindex) = 1:numel(sortindex);                           % Get the inverse of the sorting index for later use

if n_shoot == 1                                                             % Arrange integrated end values for the boundary residual
    Z_endDiagMatr = obj.perturbMatr(Z_end,dim, n_char);
else
    TempMatr = reshape(Z_end, [], dim + 1);
    Z_endDiagMatr = repmat(TempMatr(:, 1), [1, dim*n_char+1]);
end

Z0DiagMatr = obj.perturbMatr(...
    permute(Z0_perturbed(:, 1, :, :), [1, 3, 4, 2]),...                     % Select the first shooting nodes and reorder the dimensions expected by perturbMatr
    dim, n_char);                                                           % which then separates the perturbations into corresponding columns.

Z0Vec = reshape(Z0DiagMatr.', [dim*(dim*n_char + 1), n_char]);              % Reorder the data so each column contains the input values for later spline interpolation.

Z_endReshaped = reshape(Z_endDiagMatr, [dim, n_char, dim*n_char+1]);                % Reshape the endpoint sample matrix
Z_endSorted = reshape(Z_endReshaped(:, sortindex, :), [dim*n_char, dim*n_char+1]);  % and reorder it to match the end phases.

%% periodic boundary conditions
x_values = [obj.phi(2, :), 2*pi];                                    % Append 2*pi
y_values = [Z0Vec, Z0Vec(:, 1)];                                            % Append first value

spline = csape(x_values, y_values, 'periodic');

evalSpline = fnval(spline, finPhasesSorted);                                % Evaluate the interpolated polynomial of the initial states at the sorted final phases
evalSplineMatr = reshape(evalSpline, [n_char*dim + 1, n_char*dim]).';       % Reshape the vector into a matrix

ResidualMatr = Z_endSorted - evalSplineMatr;                                % Define residual of the end points of the characteristics

ResidualMatrReshaped = reshape(ResidualMatr, [dim, n_char, dim*n_char+1]);  % Reshape the Residual
ResidualMatr = reshape( ResidualMatrReshaped(:, inverseSortindex, :), ...   % Reverse the sorting
    [dim*n_char,dim*n_char+1]);

%% Assemble the complete residual matrix
if n_shoot == 1
    FullResidualMatr = ResidualMatr;                                        % For single shooting the residual matrix is defined completely by ResidualMatr.
else
    node_index = reshape(1:n_state,[dim,n_shoot,n_char]);                   % Index matrix corresponding to Z0_nodes(dim,n_shoot,n_char)
    boundary_rows = reshape(node_index(:,end,:),[],1);                      % Indices of the rows containing the quasiperiodic boundary equations

    expandedResidualMatr = obj.expandResidualMatr( ...                      % Embed the residual on the boundary into the
    ResidualMatr,dim,n_char,n_shoot);                                       % full residual matrix for the multiple shooting method.

    interior_start_nodes = reshape( ...                                     % Reshaped index vector containing the current nodes
        node_index(:,1:end-1,:),[],1);                                      % needed for the interior continuity conditions.

    interior_end_nodes = reshape( ...                                       % Indices of the following shooting nodes, which serve
        node_index(:,2:end,:),[],1);                                        % as the initial values of the subsequent intervals.

    Z0SubDiagMatr = obj.perturbMatrSub( ...                                 % Rearrange the values of Z0_perturbed in a block-diagonal structure.
        Z0_perturbed,dim,n_char,n_shoot);

    Z_SubEndDiagMatr = obj.perturbMatrSub( ...                              % Rearrange the values of Z_SubEnd in a block-diagonal structure.
        Z_SubEnd,dim,n_char,n_shoot);

    FullResidualMatr = zeros( ...                                           % Initialize the residual matrix for the multiple shooting method
        n_state,n_state+1,'like',ResidualMatr);                             % using the same data type, attributes and sparsity of ResidualMatr.

    %% Inner continuity equations
    FullResidualMatr(interior_start_nodes,:) = ...                          % Calculate the residual [Z_end(i) - Z_start(i+1) = 0]
        Z_SubEndDiagMatr(interior_start_nodes,:) ...                        % for the inner shooting nodes.
        - Z0SubDiagMatr(interior_end_nodes,:);

    %% Quasiperiodic boundary equations
    FullResidualMatr(boundary_rows,:) = ...                                 % Insert values of the residual on the boundary
        expandedResidualMatr(boundary_rows,:);                              % into the full residual matrix.

    boundaryFlowDelta = Z_SubEndDiagMatr(boundary_rows,:) ...               % Calculate the change in the integrated endpoint
        - repmat( ...                                                       % caused by perturbing the initial value of the last shooting interval.
            Z_SubEndDiagMatr(boundary_rows,1), ...
            [1,n_state+1]);

    FullResidualMatr(boundary_rows,:) = ...                                 % Complete the quasiperiodic boundary residuals
        FullResidualMatr(boundary_rows,:) + boundaryFlowDelta;              % by adding the integrated endpoint variations.
end

%% Assemble basic residual and Jacobian
f = FullResidualMatr(:,1);                                                  % Extract the unperturbed residual
f = f(:);                                                                   % and convert it into a vector.

J = sparse( (FullResidualMatr(:,2:end) ...                                  % Approximate each Jacobian column by the
            - FullResidualMatr(:,1))./dz);                                  % forward difference [R(Z + dz*e_j) - R(Z)]/dz.


%% Add parts specifically for the mixed case 
node_index = reshape(1:n_state,[dim,n_shoot,n_char]);
dP = zeros(1,n_state);
for i_shoot = 1:n_shoot
    for i_char = 1:n_char
        for i_dim = 1:dim
            deltaF = zeros(size(F_phase{i_shoot}));
            deltaF(:,i_char,:) = reshape( ...
                F_perturbed{i_shoot}(:,:,i_char,i_dim+1), ...                       % Select the correct perturbed solution
                [numel(Ik_phase{i_shoot}),1,dim]) - F_phase{i_shoot}(:,i_char,:);   % and subtract the corresponding unperturbed solution.
  
            dP(node_index(i_dim,i_shoot,i_char)) = ...                                  % Approximate the derivative of the phase condition
                obj.poincare_int(deltaF,F1_phase{i_shoot},Omega,Ik_phase{i_shoot})/dz;  % using forward difference.
        end
    end
end
J = [J,sparse(n_state,1);sparse(dP),0];

f(end+1,1) = obj.poincare_int_multiple( ...
    F_phase,F1_phase,Omega,Ik_phase);

% Differentiate the complete residual at fixed reference data, including
% phase remapping and any frequency dependence of the ODE wrapper.
domega = sqrt(eps)*(1+abs(Omega(2)));                                       % Define a small perturbation in the frequency
y_perturbed = y;
y_perturbed(end-1) = y_perturbed(end-1) + domega;                           % Only perturb the frequency in the solution
f_perturbed = obj.qp_SHM_mixed_fun(y_perturbed,DYN);                        % Calculate new residual with the perturbed frequency
J(:,end) = sparse((f_perturbed-f)./domega);                                 % Set the last column of the Jacobian
end


