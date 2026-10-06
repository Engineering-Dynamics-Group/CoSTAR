% Integral poincare phase condition for multiple shooting
% This function is a method of subclass AM_QPS_SHM.
% Calculates the integral poincare phase condition for solutions which
% require a phase condition (mixed-case, full-autonomous case)
% by summing the contributions of the individual shooting intervals.
%
%@obj:  ApproximationMethod subclass object
%@F:    Cell array of values of solution along characteristics
%@F1:   Cell array of values of derivative of reference solution along characteristics
%@Omega:Frequency vector [Omega1,Omega2]
%@Ik:   Cell array of time grids for the shooting intervals
%
%@P:    Value of integral poincare phase condition
%
% F{i} and F1{i} have size [number of time points,n_char,dim].
% Ik{i} contains the corresponding physical times, including both endpoints.
% Reference values must be provided on the same grid as the solution values
% and remain fixed during the nonlinear solution process.
% Adjacent intervals share a time endpoint, but their solution values may
% differ while the multiple shooting continuity conditions are not satisfied.
%
function P = poincare_int_multiple(obj,F,F1,Omega,Ik)

%% Initialization
n_shoot = obj.n_shoot;                                                     % Get number of shooting intervals
assert(iscell(F) && iscell(F1) && iscell(Ik), ...
    'F, F1 and Ik must be cell arrays with one entry per shooting interval.');
assert(numel(F)==n_shoot && numel(F1)==n_shoot && numel(Ik)==n_shoot, ...
    'F, F1 and Ik must contain n_shoot entries.');
P = 0;                                                                     % Initialize integral poincare phase condition

%% Integrate phase condition over all shooting intervals iteratively
for i = 1:n_shoot
    % Assert that the dimensions match
    assert(isequal(size(F{i}),size(F1{i})) ...
        && size(F{i},1)==numel(Ik{i}) ...
        && size(F{i},2)==obj.n_char && size(F{i},3)==obj.n, ...
        'Solution and reference values must match the time grid and state dimensions.');
    % Assert that the time intervals fit
    assert(isvector(Ik{i}) && numel(Ik{i})>=2 && all(diff(Ik{i})>0), ...
        'Each shooting interval must have at least two strictly increasing time points.');

    h0 = sum(F1{i}.*F{i},3);                                               % Scalar product for state-space variables (stored in the 3rd array dimension)
    P_i = trapz(obj.phi(2,:),trapz(Omega(1,1)*Ik{i},h0,1));                % Integrate phase condition numerically similar to the single shooting case
    P = P + P_i;                                                           % Add contribution of the current shooting interval
end

end
