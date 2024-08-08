clear all; close all; clc;
addpath(genpath(cd));

%% DUFFING EXAMPLE
%% Parameters and Equation
mu_limit = [0.01, 3];
epsilon = mu_limit(1);

param = {epsilon};
Fcn = @(t,z,param)vdP_auto_ap(t,z,param);

active_parameter = 1;
% non_auto_freq = @(mu) mu;
auto_freq = 1;

%% Define the input options
options.system              = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol             = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','MSHM','act_param',active_parameter);   %Properties of the solution
options.opt_cont            = costaropts('step_control','angle','mu_limit',mu_limit,'step_width',0.1,'step_width_limit',[0.001,0.3]);        %Properties for continuation
options.opt_approx_method   = costaropts('solver','ode45','n_shoot',5);
options.opt_init            = costaropts('ic',[0;1]);
options.opt_stability       = costaropts('iterate_bfp','on');

%% Continuation
[S,DYN] = costar(options);

opts = costaropts('zaxis','max2');
S.contplot(DYN,opts);



figure('Color','white')
plot(S.arclength,abs(S.multipliers),'LineWidth',2)
grid on
