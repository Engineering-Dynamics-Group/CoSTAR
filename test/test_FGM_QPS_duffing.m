%% Example: Duffing Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
D = 0.05;            ratio = sqrt(2);           % Damping factor and ratio of excitation frequencies (Omega_2/Omega_1)
f1 = 0.25;           f2 = 0.25;                 % Amplitudes of excitation
kappa = 1.5;                                    % Coefficient of nonlinear stiffness

mu_limit = [0.8, 2];   eta = mu_limit(2);       % Limits of continuation diagram and mu-value at start of continuation

param = {D,kappa,f1,f2,eta,ratio};
active_parameter = 5;

Fcn = @(t,z,param)duffing_ap_qp(t,z,param);
non_auto_freq = @(mu) [mu,ratio.*mu];

% Init data
K3 = [0, 1, 0; 0, 0, 1];
% K3 = [0 1 0 3 0 -1 1 -2 2 -3  -2   2 -4  -5  -4  -5  3  3  4  5  7;...        %Higher harmonics will be guessed
%       0 0 1 0 3  2 2  1 1  2   3   3  3   2   1   6  2  4  1  0  0];

% c_max =     [ 0.1188,    0.1443;          % Fourier coefficients for mu = 0.8
%              -0.0104,   -0.0075];
% s_max =     [-0.1049,   -0.0539;
%              -0.0119,   -0.0205];
%
c_max2 = [  0,     0;                       % Fourier coefficients for mu = 2
          -0.2, -0.07];
s_max2 = [-0.1, -0.025;
            0,     0  ];
% c_max2 = [c_max, 0.1.*ones(2,14)];
% s_max2 = [s_max, 0.1.*ones(2,14)];


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                   % Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','quasiperiodic',...
                              'approx_method','fourier-galerkin','act_param',active_parameter,'display','error-control');   % Properties of the solution
options.opt_init = costaropts('C0',zeros(2,1),'cmatrix',c_max2,'smatrix',s_max2,'hmatrix',K3);
options.opt_approx_method = costaropts('n_FFT',2^3,'error_limit',[1e-3 1e-2]);                                              % Properties for approx_method
options.opt_cont = costaropts('pred','parable','mu_limit',mu_limit,'direction',-1);                                         % Properties for continuation
options.opt_stability = costaropts('iterate_bfp','on','n_char_st',50,'n_map',5e3);                                          % Properties for stability


%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
