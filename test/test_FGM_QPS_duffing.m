%% Example: Duffing Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
mu_limit = [0.8,2];

D = 0.05;
kappa = 1.5;
f1 = 0.25;
f2 = 0.25;
eta = mu_limit(2);
ratio = sqrt(2);

param = {D,kappa,f1,f2,eta,ratio};
active_parameter = 5;

Fcn = @(t,z,param)duffing_ap_qp(t,z,param);
non_auto_freq = @(mu) [mu,ratio.*mu];


% Init data
K3 = [0, 1, 0; 0, 0, 1];

c_max2 = [  0,     0;                      % Fourier coefficients for mu = 2
          -0.2, -0.07];
s_max2 = [-0.1, -0.025;
            0,     0  ];

% c_max =     [ 0.1188,    0.1443;
%              -0.0104,   -0.0075];
% s_max =     [-0.1049,   -0.0539;
%              -0.0119,   -0.0205];
%
% K3 = [0 1 0 3 0 -1 1 -2 2 -3  -2   2 -4  -5  -4  -5  3  3  4  5  7;...        %Higher harmonics will be guessed
%       0 0 1 0 3  2 2  1 1  2   3   3  3   2   1   6  2  4  1  0  0];
%
% c_max2 = [c_max, 0.1.*ones(2,14)];
% s_max2 = [s_max, 0.1.*ones(2,14)];


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                   % Properties of the System
options.opt_sol  = costaropts('stability','off','cont','on','non_auto_freq',non_auto_freq,'sol_type','quasiperiodic',...
                              'approx_method','fourier-galerkin','act_param',active_parameter,'display','step-control');    % Properties of the solution
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',c_max2,'Smatrix',s_max2,'Hmatrix',K3);
options.opt_approx_method = costaropts('n_FFT',2^6,'error_limit',[1e-3 1e-2]);                                              % Properties for approx_method
options.opt_cont = costaropts('pred','parable','mu_limit',mu_limit,'direction',-1);                                         % Properties for continuation

timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
