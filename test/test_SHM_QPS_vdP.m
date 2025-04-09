%% Example: van der Pol Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
epsilon = 0.3;
s = 0.4;

IC = 0.1.*ones(2,1);                                                                         %initial condition
mu_limit = [1.3,1.45];
non_auto_freq = @(mu) mu;                                                                 %non-autonomous frequencies
auto_freq = 1;

param = {epsilon,s,mu_limit(1,2)};                                            %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 3;                                                                   %Index of active parameter
Fcn = @(t,z,param)vdP_qp(t,z,param);                                        %Right-hand-side of ODE

% load('workspace_test_vdp_QPS_shooting.mat')       % OLD. Throws warning 'Could not find appropriate function on path loading function handle C:\Users\Admin\Desktop\FG-Code_aktuell\v2.1.1.15\test\test_Shooting_vdP_QPS.m>@(mu)mu' after test files have been renamed, altough path was not existing even before
load('workspace_test_SHM_QPS_vdp.mat')              %only contains s0 since DYN_init stored in workspace_test_vdp_QPS_shooting throws warning

% C1_mat = [-0.5, 0;                              % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
%           0, -1.5];
% S1_mat = [0, -1.5;
%           0.6, 0];


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'auto_freq',auto_freq,'sol_type','quasiperiodic','approx_method','shooting','act_param',active_parameter); %Properties of the solution
options.opt_cont = costaropts('step_control','corrector_iterations','pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit);                                           %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_char',50);                                                                             %Properties for sol_method (e.g. Shoot)
options.opt_init = costaropts('iv',s0); 
% options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);
options.opt_stability = costaropts('iterate_bfp','on','n_char_st',50,'n_map',1e4);                                      %Properties for stability

options.opt_cont.step_width = 0.5;                                                        %Initial step width
options.opt_cont.direction = -1;                                                        %Changes the direction of continuation 


%% Continuation
[S,DYN] = costar(options);  


%% Test Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
