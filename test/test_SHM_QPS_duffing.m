%% Example: Duffing Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
% "Difficult" parameters
% D = 0.1;              ratio = 1/sqrt(2);              % Damping factor and ratio of excitation frequencies (Omega_2/Omega_1)
% s1 = 5;               s2 = 5;                         % Amplitudes of excitation
% kappa = 0.2;                                          % Coefficient of nonlinear stiffness
% mu_limit = [3,5];                                     % Limits of continuation diagram

% "Easy" parameters - advantageous for testing stability computation
D = 0.05;            ratio = sqrt(2);                   % Damping factor and ratio of excitation frequencies (Omega_2/Omega_1)
s1 = 0.25;           s2 = 0.25;                         % Amplitudes of excitation
kappa = 1.5;                                            % Coefficient of nonlinear stiffness
mu_limit = [0.8, 2];                                    % Limits of continuation diagram

param = {D,kappa,s1,s2,mu_limit(2),ratio};              %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 5;                                   %Index of active parameter

Fcn = @(t,z,param)duffing_qp(t,z,param);                %Right-hand-side of ODE
non_auto_freq = @(mu) [mu mu*ratio];                    %non-autonomous frequencies

IC = 0.1.*ones(2,1);                                    %initial condition

% load('workspace_test_duffing_QPS_shooting.mat')       % OLD. Throws warning 'Could not find appropriate function on path loading function handle C:\Users\Admin\Desktop\FG-Code_aktuell\v2.1.1.15\test\test_Shooting_duffing_QPS.m>@(mu)[mu,mu*Omega]' after test files have been renamed, altough path was not existing even before
load('workspace_test_SHM_QPS_duffing.mat')              %Contains s0 of difficult and easy parameters at mu_limit(2) since DYN_init stored in workspace_test_duffing_QPS_shooting throws warning

% C1_mat = [  0,     0,   0;                      % Fourier coefficients to create an initial value for fsolve to find the first point on the curve at mu = 2
%           -0.2, -0.07, 0];
% S1_mat = [-0.1, -0.025, 0;
%             0,     0,   0];

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                               %Properties of the system
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','quasiperiodic','approx_method','shooting','act_param',active_parameter,'display','iter-detailed'); %Properties of the solution
options.opt_cont = costaropts('step_control','angle','pred','cubic','subspace','pseudo-arc','mu_limit',mu_limit);       %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_char',50);                                                   %Properties for the approximation method
options.opt_init = costaropts('iv',s_easy_mu_2);                                                                        %Properties for initial value
% options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);
options.opt_stability = costaropts('iterate_bfp','on','n_char_st',75,'n_map',1e3);                                      %Properties for stability

options.opt_cont.step_width = 0.25;                     %Initial step width    
options.opt_cont.step_control_param = [2, 5];           %Step control parameters
options.opt_cont.direction = -1;                        %Changes the direction of continuation 


%% Continuation
timer = tic;                                            % Record current time
[S,DYN] = costar(options);                              % Calculate initial solution and continue the curve
time = toc(timer);                                      % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
