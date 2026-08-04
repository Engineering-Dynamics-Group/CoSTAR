% This Continuation method calculates the initial slope, which is used for secant, parble or ... 
% cubic predictor (if they were chosen) to determine the direction vector if only 1 curve point exists.
% 
% @obj:     continuation class object

function obj = initial_slope(obj,DYN,AM)

y0 = obj.y0;                                                            % Initial solution
de = 1e-4;                                                              % Small value to pertub the mu-value of y0
y0_de = [y0(1:end-1); y0(end)+obj.direction*de];                        % New y0 value to compute a second solution

AM.IF_up_res_data(y0,DYN);                                              % Update AM properties and set y0 as initial value
Fcn = @(y) AM.res_fun_init(y,y0_de);                                    % Function wrapper to set the complete residuum function

[ys,~,secant_flag,~,~] = fsolve(Fcn,y0,obj.fsolve_opts);                % Compute solution

if secant_flag < 1                                                      % No solution found
    obj.p_use_tangent = true;                                           % Use tangent to get the first direction vector
    info_text = 'No intermediate curve point found! Using tangent as direction vector for first step.';
else
    obj.p_initial_slope = (ys-y0);                                      % Compute initial slope
    info_text = 'Intermediate curve point found! Using secant as direction vector.';
end

disp(info_text);
write_log(DYN,info_text);
obj.p_last_msg = sprintf('%s\n',info_text);

end