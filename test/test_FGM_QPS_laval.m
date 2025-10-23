%%    Example: Laval rotor (quasi-periodic)     %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
load('workspace_init_laval_QPS_FGM');

Di = 0.2;     Delta = 1/3;     d3 = 0.25;       % Damping parameters
e = 0.25;                                       % Eccentricity
Fg = 0.3924;                                    % Weight force

mu_limit = [1.72,2.5];    eta = mu_limit(2);    % Limits of continuation diagram and mu-value at start of continuation | NSB at mu = 1.71477

param = {eta, Di, Delta, e, d3, Fg};            % Parameter array
active_parameter = 1;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) mu;                       % Non-autonomous frequency
Fcn = @(t,z,param) laval_qp(t,z,param);         % Right hand side of ODE


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);    %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'non_auto_freq',non_auto_freq,'sol_type','qps','approx_method','fgm','act_param',active_parameter,'display','error-control');      %Properties of the solution
options.opt_cont = costaropts('step_control','angle','direction',-1,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.02,'max_cont_step',1e4);                                                             %Properties for continuation
options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','int_poincare');   %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('c0',zeros(4,1),'cmatrix',c_max,'smatrix',s_max,'hmatrix',K3);
options.opt_stability = costaropts('iterate_bfp','on','n_char_st',75,'n_map',1e4);                                      %Properties for stability


%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
