%% Example: Coupled van der Pol Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
epsilon =   0.1;
alpha = 0.1;
beta = 1.1;
% param = {eps,alpha,beta};

IC = 1.*ones(4,1);                                                                         %initial condition
mu_limit = [0.1,0.4]; 
auto_freq = [1.04,1.49];            

mu_start = mu_limit(1,1);

param = {mu_start,alpha,beta};                                            %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 1;                                                     %Index of active parameter
Fcn = @(t,z,param)coupledvdp(t,z,param);                                  %Right-hand-side of ODE

% load('workspace_test_cvdP_QPS_shooting.mat')       % OLD. Throws warning 'Could not find appropriate function on path loading function handle C:\Users\Admin\Desktop\FG-Code_aktuell\v2.1.1.15\test\test_Shooting_cvdP_QPS.m>@(t,z,param)coupledvdp(t,z,param)' after test files have been renamed, altough path was not existing even before
load('workspace_test_SHM_QPS_cvdP.mat')              %only contains s0 since DYN_init stored in workspace_test_cvdp_QPS_shooting throws warning


% C1_mat = [0,   0,   0;                          % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
%           0, -1.4,  0;
%           2,   0,   0;
%           0,  2.1,  0];
% S1_mat = [2,   0,   0;
%           0,  1.4,  0;
%           0,   0,   0;
%           0,  2.1,  0];

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);                                                                          %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','quasiperiodic','approx_method','shooting','act_param',active_parameter,'display','full'); %Properties of the solution
options.opt_cont = costaropts('step_control','off','pred','parable','subspace','pseudo-arc','mu_limit',mu_limit);              %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_char',35);                                                                             %Properties for sol_method (e.g. Shoot)
options.opt_init = costaropts('iv',s0); 
% options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);

options.opt_cont.step_width = 0.1;                                                       %Initial step width                                                                                                                                                                                                                     
options.opt_cont.direction = 1;                                                          %Changes the direction of continuation 


%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Test Postprocessing
% benchmark_postprocess_quasiperiodic(DYN,S);
