%%    Example: Laval rotor (quasi-periodic)     %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
Di = 0.2;     Delta = 1/3;     d3 = 0.25;       % Damping parameters
e = 0.25;                                       % Excentricity
Fg = 0.3924;                                    % Weight force

mu_limit = [1.72,2.5];    mu0 = mu_limit(2);    % Limits of continuation diagram and mu-value at start of continuation | NSB at mu = 1.71477

param = {mu0, Di, Delta, e, d3, Fg};            % Parameter array
active_parameter = 1;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) mu;                       % Non-autonomous frequency
Fcn = @(t,z,param) laval_qp(t,z,param);         % Right hand side of ODE

auto_freq = 1;                                  % Initial value of autonomous frequency
c0 = [-0.2; -0.2; 0; 0];                        % Fourier coefficient used to create an initial value for fsolve    
C1_mat = [  0,   1.35, 0;                       % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
          -0.35,  0,   0;
          -0.95,  0,   0;
            0,  -1.35, 0];
S1_mat = [-0.38,   0,   0;
            0,   -1.35, 0  
            0,   -1.35, 0
          0.875,   0,   0];


%% Properties
options.system   = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param);                                                           % Properties of the system
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','off',...      % Properties of the solution
                              'auto_freq',auto_freq,'non_auto_freq',non_auto_freq,'act_param',active_parameter);                    % Properties of the solution
options.opt_init = costaropts('c0',c0,'c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                       % Properties for initial solution
options.opt_approx_method = costaropts('n_int_1',30,'scheme_1','central','approx_order_1',6,...                                     % Properties of approximation method FDM
                                       'n_int_2',30,'scheme_2','central','approx_order_2',6);                                       % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','secant','direction',-1,'display','step_control_info');                    % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.5;
options.opt_cont.step_control = 'angle';
%options.opt_cont.step_control_param = [3, 5/180*pi]; 


%% Continuation
tic                                             % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
toc                                             % Display elapsed time since tic 


%% Comparison with FGM
%{
load('workspace_init_laval_QPS_FGM');

options.system   = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param);                                                           % Properties of the system
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off',...       % Properties of the solution
                              'auto_freq',auto_freq,'non_auto_freq',non_auto_freq,'act_param',active_parameter);                    % Properties of the solution
options.opt_init = costaropts('c0',zeros(4,1),'cmatrix',c_max,'smatrix',s_max,'hmatrix',K3);                                        % Properties for initial solution
options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','int_poincare');                                                     % Properties of approximation method FGM
options.opt_cont = costaropts('mu_limit',mu_limit,'maxcontstep',1e4,'step_width',0.01,'direction',-1);                              % Properties for continuation

tic
[S1,DYN1] = costar(options);
toc                                                                   % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Plot continuation 1 (using FGM) into figure of continuation plot FDM

% benchmark_postprocess_quasiperiodic(DYN1,S1)
%}


%% Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);     % Postprocessing
