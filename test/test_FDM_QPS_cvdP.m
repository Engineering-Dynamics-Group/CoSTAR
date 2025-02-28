%%  Example: Coupled van der Pol Oscillator (quasi-periodic)  %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
alpha = 0.1;    beta = 1.1;

mu_limit = [0.1, 1];                          % Limits of continuation diagram (mu > 1.3 and fsolve calculates J using central FD: no convergence sometimes ?)
epsilon0 = mu_limit(1);                         % mu-value at start of continuation

param = {epsilon0, alpha, beta};                % Parameter array
active_parameter = 1;                           % Defines where the continuation parameter is located within param

Fcn = @(t,z,param) coupledvdp(t,z,param);

auto_freq = [1.04, 1.5];                        % Initial autonomous frequencies
C1_mat = [0,   0,   0;                          % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
          0, -1.4,  0;
          2,   0,   0;
          0,  2.1,  0];
S1_mat = [2,   0,   0;
          0,  1.4,  0;
          0,   0,   0;
          0,  2.1,  0];


%% Properties
options.system   = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param,'info','continuation of coupled van der Pol');              % Properties of the System
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','on',...       % Properties of the solution
                              'auto_freq',auto_freq,'act_param',active_parameter,'display','step-control');                         % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                               % Properties for initial solution
options.opt_approx_method = costaropts('n_int_1',35,'points_1',[-4,-3,-2,-1,0,1,2],...                                              % Properties of approximation method FDM
                                       'n_int_2',35,'points_2',[-4,-3,-2,-1,0,1,2]);                                                % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','secant');                                                                 % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 1.5;
options.opt_cont.step_control = 'off';
% options.opt_cont.step_control_param = [2, 5]; 


%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Comparison with FGM

%{
load('workspace_init_cvdp_QPS_FGM');
    
K3 = [0 1 0  3  0 -1 1;...      % Readjust the higher harmonics matrix
      0 0 1  0  3  2 2];

options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);                                                           % Properties of the System
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off',...       % Properties of the solution
                              'auto_freq',auto_freq,'act_param',active_parameter);                                                  % Properties of the solution
options.opt_init = costaropts('c0',zeros(4,1),'cmatrix',c_max,'smatrix',s_max,'hmatrix',K3);                                        % Properties for initial solution
options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','int_poincare','n_hh_max',50,'error_limit',[1e-3,1e-2]);             % Properties for approx_method
options.opt_cont = costaropts('mu_limit',mu_limit,'direction',1,'step_width',0.05,'maxcontstep',1e4);                               % Properties for continuation

tic                                                                     % Record current time
[S1,DYN1] = costar(options);                                            % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Plot continuation 1 (using FGM) into figure of continuation plot FDM

benchmark_postprocess_quasiperiodic(DYN1,S1)
%}


%% Postprocessing
benchmark_postprocess_quasiperiodic(DYN,S);
