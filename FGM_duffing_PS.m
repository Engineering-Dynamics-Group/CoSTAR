clear all; close all; clc;
addpath(genpath(cd));

%% Parameters
kappa = 0.3;
D = 0.05;
g = 1;

%% DUFFING EXAMPLE
IC = [1.5;0];     mu_limit = [0.1, 2.5];
non_auto_freq = @(mu) mu;

eta = mu_limit(1);     param = {kappa,D,eta,g};     active_parameter = 3;
Fcn = @(t,z,param)duffing_ap(t,z,param);

%% Continuation
load('workspace_test_duffing_PS_FGM');
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','fourier-galerkin','act_param',active_parameter);   %Properties of the solution
options.opt_cont = costaropts('step_control','angle','direction',1,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.1);                 %Properties for continuation
options.opt_approx_method = costaropts('n_FFT',2^6,'error_control','on','error_limit',[0,0.001]);                                                                                                           %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',[c_max,[0.01;0.01]],'Smatrix',[s_max,[0.01;0.01]],'Hmatrix',[0,1,3]);
options.opt_stability       = costaropts('iterate_bfp','on','solver','ode45','n_shoot',7);

[S,DYN] = costar(options);

figure('Color','white')
plot(S.arclength,abs(S.multipliers))