% Time grids for the integral poincare phase condition
% This function is a method of subclass AM_QPS_SHM.
% Generates time grids for the individual shooting intervals and evaluates
% the fixed reference gradient at corresponding normalized time positions.
%
%@obj:      ApproximationMethod subclass object
%@T:        Current integration period
%@Ik_phase: Cell array of physical time grids
%@tau_phase:Cell array of normalized time grids
%@F1_phase: Cell array of reference gradients (optional output)
%
function [Ik_phase,tau_phase,F1_phase] = getPhaseGrids(obj,T)

tau = linspace(0,1,obj.reso_phase);                                        % Normalized time grid for phase quadrature
Ik_phase = cell(1,obj.n_shoot);
tau_phase = cell(1,obj.n_shoot);
F1_phase = cell(1,obj.n_shoot);

for i = 1:obj.n_shoot
    bounds = [(i-1)/obj.n_shoot,i/obj.n_shoot];                             % Set normalized start and end point of the i-th shooting interval 
    tau_i = [bounds(1),tau(tau>bounds(1) & tau<bounds(2)),bounds(2)];       % Array of all tau values that lie inside the closed interval from bound_1 to bound_2
    
    if numel(tau_i)==2
        tau_i = [bounds(1),mean(bounds),bounds(2)];                         % Add a midpoint so the ODE solver returns values on the requested grid
    end                                                                     

    tau_phase{i} = tau_i;                                                   % Save the normalized times for the evaluation
    Ik_phase{i} = T*tau_i;                                                  % Also save the positions in the actual time t

    % If only 2 outputs are requested, the selection and interpolation 
    % of the gradient is skipped
    if nargout < 3
        continue
    end

    %Calculation of the gradient
    if isempty(obj.phase_reference)                                         % Check if new phase_reference structure is present.
        F1 = obj.Y_old{1,3};                                                
        tau_ref = linspace(0,1,size(F1,1));                                 % If that isn't the case, normalized times are assigned to the stored reference gradient samples (F1)
    else
        ref = obj.phase_reference;
        if isscalar(ref.tau)
            source = 1;                                                     % If only a single interval in hypertime is given, use it for every shooting interval
        else
            assert(numel(ref.tau)==obj.n_shoot,'Rebuild the phase reference after changing n_shoot.');  %If multiple intervals are given, the amount has to match n_shoot
            source = i;                                                     % Use i-th interval for the gradient selection
        end

        F1 = ref.gradient{source};                                          % Select the gradient corresponding to the chosen interval
        tau_ref = ref.tau{source};                                          % Select the normalized time grid of the chosen reference segment
    end

    F1_mat = reshape(F1,size(F1,1),[]);                                     % Reshape the gradients of the reference solution for later interpolation
    F1_interpolated = interp1(tau_ref,F1_mat,tau_i,'linear');               % Interpolate every column from tau_ref to tau_i
    F1_phase{i} = reshape(F1_interpolated,[numel(tau_i),obj.n_char,obj.n]); % Restore the dimensions that represent [evaluation times, characteristics, state components]
end
end
