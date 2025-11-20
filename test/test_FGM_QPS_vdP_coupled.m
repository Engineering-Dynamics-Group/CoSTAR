%%  Example: Coupled van der Pol Oscillator (quasi-periodic)  %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
eps_limit = [0.01,1];

alpha = 0.1;
beta = 1.1;
param = {eps_limit(1),alpha,beta};

cont = 1;
active_parameter = 1;
Fcn = @(t,z,param) vdP_coupled(t,z,param);

load('workspace_init_cvdp_QPS_FGM');

%Readjust the higher harmonics matrix
K3 = [0 1 0  3  0 -1 1;...
    0 0 1  0  3  2 2];


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);    %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','qps','approx_method','fgm','act_param',active_parameter,'display','error-control');      %Properties of the solution
options.opt_cont = costaropts('step_control','angle','direction',1,'pred','tangent','subspace','pseudo-arc','mu_limit',eps_limit,'step_width',0.05,'max_cont_step',1e4);               %Properties for continuation
options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','int_poincare','n_hh_max',50);   %Properties for approx_method (e.g. Shoot)
options.opt_init = costaropts('c0',zeros(4,1),'cmatrix',c_max,'smatrix',s_max,'hmatrix',K3);
options.opt_stability = costaropts();                                      %Properties for stability



%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
