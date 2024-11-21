% This Continuation method calculates the initial slope, which is used for secant, parble or ... 
% cubic predictor (if they were chosen) to determine the direction vector if only 1 curve point exists.
% 
% @obj:     continuation class object

function obj = initial_slope(obj,DYN,AM)

newtonOpts = optimoptions('fsolve','Display','iter-detailed','MaxFunEvals',1e5,'FiniteDifferenceType','forward');
y0 = obj.y0;
de = 1e-4;

if(strcmpi(DYN.sol_type,'quasiperiodic')&&strcmpi(DYN.approx_method,'shooting'))
    Fcn = @(y)AM.fun_Jac_wrapper_init(y,y0(end)+obj.direction.*de,DYN);                        % Function wrapper for initial solution, if Jacobian is supplied
    newtonOpts.SpecifyObjectiveGradient = true;
elseif strcmpi(DYN.approx_method,'finite-difference')                   % Special corrector function due to specification of Jacobian matrix
    newtonOpts.SpecifyObjectiveGradient = true;                         % newtonOpts.CheckGradients = true; can be used to automatically check the Jacobian matrix -> Since R2023b: checkGradients is recommended
    Fcn = @(y) AM.corr_fun_init_FDM(y,y0(end)+obj.direction.*de);
else
    Fcn = @(y)[AM.res(y);y(end)-(y0(end)+obj.direction.*de)];                               % Last entry is not necessary, but this way the Jacobian has the correct dimension
end

[ys,~,secant_flag,~,~] = fsolve(Fcn,y0,newtonOpts);
obj.p_initial_slope = (ys-y0);

if(secant_flag~=1)
    obj.p_use_qr = true;
    info_text = 'No second curve point found! Switched to tangent predictor for initial step!';
else
    info_text = 'Second cuve point found!';
end

disp(info_text); disp(' ');
write_log(DYN,info_text);

end