% Function sets new calculated curve point as active curve point to calculate next curve point
%
% @obj: Continuation class object

function obj = iterate_data(obj)

    obj.p_arcl_0        = obj.p_arcl_1;                     % Arclength -> This must be called before obj.y0 = obj.p_y1 !!!
    obj.p_mu0           = obj.p_y1(end,1);                  % Bifurcation parameter
    obj.p_y0_old        = [obj.p_y0_old(2:end), num2cell(obj.y0,1)];  % The 3 last "old" curve points, stored in cell arrays
    obj.y0              = obj.p_y1;                         % Solution vector
    obj.p_dy_old        = obj.dy0;                          % Direction vector
    obj.p_J0            = obj.p_J1;                         % Jacobian matrix
    obj.p_it            = obj.p_output.iterations;          % Number of corrector iterations for step control
    obj.p_r_old         = obj.p_r;                          % Factor which adapts step width
    obj.p_n_unstable_0  = obj.p_n_unstable_1;               % Number of unstable multipliers
  
    if strcmpi(obj.step_control,'pid')
        obj.p_e_old_old  = obj.p_e_old;                     % Used for PID step control
        obj.p_e_old      = obj.p_e;                         % Used for PID step control
    end
    
end