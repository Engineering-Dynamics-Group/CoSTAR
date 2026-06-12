% This Continuation method calculates the initial slope, which is used for secant, parble or ... 
% cubic predictor (if they were chosen) to determine the direction vector if only 1 curve point exists.
% 
% @obj:     continuation class object

function obj = initial_slope(obj,DYN,AM)

y0 = obj.y0;
de = 1e-4;

if strcmpi(DYN.approx_method,'shooting')
    obj.fsolve_opts.MaxIter = 50;
    AM.IF_up_res_data(obj,DYN);                                         % Pass information to the ApproxMethod object
    Fcn = @(y) AM.fun_Jac_wrapper_init(y,y0(end)+obj.direction.*de);    % Function wrapper for initial solution, if Jacobian is supplied
    obj.fsolve_opts.SpecifyObjectiveGradient = true;                    % Jacobian matrix is passed by the user
elseif strcmpi(DYN.approx_method,'finite-difference')                   % Special corrector function due to specification of Jacobian matrix
    AM.IF_up_res_data(obj);                                             % Pass information to the ApproxMethod object
    obj.fsolve_opts.SpecifyObjectiveGradient = true;                    % Jacobian matrix is passed by the user
    Fcn = @(y) AM.corr_fun_init_FDM(y,y0(end)+obj.direction.*de);
else
    AM.IF_up_res_data(obj);                                             % Pass information to the ApproxMethod object
    Fcn = @(y)[AM.res(y);y(end)-(y0(end)+obj.direction.*de)];           % Last entry is not necessary, but this way the Jacobian has the correct dimension
end

[ys,~,secant_flag,~,~] = fsolve(Fcn,y0,obj.fsolve_opts);
obj.p_initial_slope = (ys-y0);

if secant_flag < 1                                                      % No solution found
    obj.p_use_tangent = true;                                           % Use tangent to get the first direction vector
    info_text = 'No intermediate curve point found! Using tangent as direction vector for first step.';
else
    info_text = 'Intermediate curve point found! Using secant as direction vector.';
end

disp(info_text);
write_log(DYN,info_text);
obj.p_last_msg = sprintf('%s\n',info_text);

end