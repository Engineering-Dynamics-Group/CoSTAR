% Function checks whether calculated curve point is in the boundaries given
% by mu_limit and the maximal number of iterations is not surpassed. It
% furthermore dispays informations in the console
%
%@obj:      Continuation class object
function obj = check_limits(obj)

obj.p_local_cont_counter = obj.p_local_cont_counter + 1;    % Set counter to next integer

formatSpec = '%.4f';

% Display information about curve point
if strcmpi(obj.display,'on') || strcmpi(obj.display,'step_control_info')
    disp(append('Iter: ',num2str(obj.p_local_cont_counter),' -- mu = ',num2str(obj.p_mu0,formatSpec),' -- stepwidth = ',num2str(obj.step_width,formatSpec)));
end

% Check the stopping criteria
if(obj.p_mu0>obj.mu_limit(1,2))         % Check if calculated curve point is above upper limit of mu
    obj.p_contDo = 0;
    clear fsolve;
    disp(['Maximal value of continuation parameter mu_max = ', num2str(obj.mu_limit(1,2)), ' reached: mu = ', num2str(obj.p_mu0)]);
    obj.p_stopping_flag = 1;
elseif(obj.p_mu0<obj.mu_limit(1,1))     % Check if calculated curve point is below bottom limit of mu
    obj.p_contDo = 0;
    clear fsolve;
    disp(['Minimal value of continuation parameter mu_min = ', num2str(obj.mu_limit(1,1)), ' reached: mu = ', num2str(obj.p_mu0)]);
    obj.p_stopping_flag = 1;
elseif(obj.p_local_cont_counter>=obj.max_cont_step)      % Check if maximal number of continuation points has been surpassed
    obj.p_contDo = 0;
    clear fsolve;
    disp(['Maximal number of continuation steps max_cont_step = ', num2str(obj.max_cont_step), ' reached.']);
    obj.p_stopping_flag = 2;
end

end