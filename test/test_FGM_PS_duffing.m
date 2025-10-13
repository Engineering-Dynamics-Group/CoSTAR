%%    Example: Duffing Oscillator (periodic)   %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
c = 1;
kappa = 0.3;
D = 0.05;
g = 1;

IC = [1.5;0];     mu_limit = [0.01,2.5];
non_auto_freq = @(mu) mu;

eta = mu_limit(1);     param = {kappa,D,eta,g,c};     active_parameter = 3;
Fcn = @(t,z,param)duffing_ap(t,z,param);


%% Single Solution
load('workspace_test_duffing_PS_FGM');
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','off','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','fourier-galerkin','act_param',active_parameter);   %Properties of the solution
options.opt_approx_method = costaropts('n_FFT',2^6);                                                                                                           %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',[c_max,[0.01;0.01]],'Smatrix',[s_max,[0.01;0.01]],'Hmatrix',[0,1,3]);
options.opt_stability       = costaropts('iterate_bfp','on');

[S,DYN] = costar(options);                                                                                                                              %Calculate initial solution and continue the curve to set limits


%% Continuation 1
load('workspace_test_duffing_PS_FGM');
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol  = costaropts('display','iter-detailed','stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','fourier-galerkin','act_param',active_parameter);   %Properties of the solution
options.opt_cont = costaropts('step_control','angle','direction',1,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.1);                 %Properties for continuation
options.opt_approx_method = costaropts('n_FFT',2^6,'error_control','off');                                                                                                           %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',c_max,'Smatrix',s_max,'Hmatrix',[0,1,3]);
options.opt_stability       = costaropts('iterate_bfp','on');

[S1,DYN1] = costar(options);                                                                                                                              %Calculate initial solution and continue the curve to set limits


%% Continuation 2
load('workspace_test_duffing_PS_FGM');
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);    %Properties of the System
options.opt_sol  = costaropts('display','iter','stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','fourier-galerkin','act_param',active_parameter);   %Properties of the solution
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.05,'step_control','angle','direction',1,'pred','tangent','subspace','pseudo-arc');                 %Properties for continuation
options.opt_approx_method = costaropts('n_FFT',2^8,'error_control','on','error_limit',[1e-4,1e-3]);                                                                                                           %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',[c_max,[0.01;0.01]],'Smatrix',[s_max,[0.01;0.01]],'Hmatrix',[0,1,3]);
options.opt_stability       = costaropts('iterate_bfp','on');

timer = tic;                                    % Record current time
[S2,DYN2] = costar(options);                    % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic
% open(append('CoSTAR_Log_',DYN2.DYN_id,'.txt'))


%% Test Postprocessing
benchmark_postprocess_periodic(DYN1,S1);
benchmark_postprocess_periodic(DYN2,S2);
