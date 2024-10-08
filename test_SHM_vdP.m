clear all; close all; clc;
addpath(genpath(cd));

%% Van-der-Pol EXAMPLE
%% Parameters and Equation
mu_limit = [0.5, 3];
epsilon = 0.2;
s = 0.2;
eta = mu_limit(2);

param = {epsilon,s,eta};
Fcn = @(t,z,param)vdP_qp(t,z,param);

active_parameter = 3;
non_auto_freq = @(mu) mu;
dir = -1;


%% Define the input options
options.system              = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol             = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','SHM','act_param',active_parameter);   %Properties of the solution
options.opt_cont            = costaropts('step_control','angle','mu_limit',mu_limit,'step_width',0.1,'step_width_limit',[0.001,0.3],'direction',dir);        %Properties for continuation
options.opt_approx_method   = costaropts('solver','ode45','n_shoot',6);
options.opt_init            = costaropts('ic',[0;1]);
options.opt_stability       = costaropts('iterate_bfp','on');

%% Continuation
[S,DYN] = costar(options);

opts = costaropts('zaxis','max2');
S.contplot(DYN,opts);

opt_solget = costaropts('space','trajectory','zaxis',@(z) z(:,1),'xaxis',@(z) z(:,2),'index',1);
S.solplot(DYN,opt_solget);



% 
% figure('Color','white')
% plot(S.arclength,abs(S.multipliers))

