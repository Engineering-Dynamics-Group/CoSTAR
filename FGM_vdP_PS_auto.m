clear all; close all; clc;
addpath(genpath(cd));

mu_limit = [0.01,3];
auto_freq = 1;   %Period-function for constant period

epsilon = mu_limit(1);     param = {epsilon};     active_parameter = 1;
Fcn = @(t,z,param)vdP_auto_ap(t,z,param);

c_max = [2    ,0.1, 0.05,0.01,0.001,0.001,0.001;     0,-0.15,-0.05,-0.01,-0.001,-0.001,-0.001];
s_max = [-1.98,0.11,0.05,0.01,0.001,0.001,0.001;    -0.007,-0.05,-0.05,-0.01,-0.001,-0.001,-0.001];


%% Continuation
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','fourier-galerkin','act_param',active_parameter);       %Properties of the solution
options.opt_cont = costaropts('step_control','angle','pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.1);             %Properties for continuation
options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','poincare');                                                                                                       %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('cmatrix',c_max,'smatrix',s_max,'c0',zeros(2,1),'Hmatrix',[0,1,3,5,7,9,11,13,15,17]);
options.opt_stability       = costaropts('iterate_bfp','on');

[S,DYN] = costar(options);                                                                                                      %Calculate initial solution and continue the curve to set limits


figure('Color','white')
plot(S.arclength,abs(S.multipliers),'LineWidth',2)
grid on
