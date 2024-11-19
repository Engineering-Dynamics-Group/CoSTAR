%% Example: Forced van der Pol Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Continuation of excitation amplitude
epsilon = 0.5;                                  % Non-linear damping coefficient
Omega = 1.2;                                    % Excitation frequency | Omega = 1.1: Two FB,the "upper" FB being at mu = 0.453789 | Omega = 1.2: NSB at mu = 0.761537

mu_limit = [0.01,0.5];    mu0 = mu_limit(1);    % Limits of continuation diagram and mu-value at start of continuation. ATTENTION: Be aware of the bifurcations!

param = {epsilon, mu0, Omega};                  % Parameter array
active_parameter = 2;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) Omega;                    % Non-autonomous frequency
Fcn = @(t,z,param) vdP_qp(t,z,param);           % Right hand side of ODE

auto_freq = 1;                                  % Initial value of autonomous frequency
C1_mat = [0, 2, 0;                              % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
          0, 0, 0];
S1_mat = [0, 0, 0;
          0, -2, 0];


% Properties
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param);                                                       % Properties of the system
options.opt_sol  = costaropts('sol_type','qps','approx_method','fdm','cont','on','stability','off','auto_freq',auto_freq,...    % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter,'display','step-control');             % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                           % Properties for initial solution
options.opt_approx_method = costaropts('n_int_1',50,'points_1',[-4,-3,-2,-1,0,1,2],...                                          % Properties of approximation method FDM
                                       'n_int_2',50,'points_2',[-4,-3,-2,-1,0,1,2]);                                            % Properties of approximation method FDM
% options.opt_approx_method = costaropts('n_int_1',50,'scheme_1','central','approx_order_1',6,...                               % Properties of approximation method FDM
%                                        'n_int_2',50,'scheme_2','central','approx_order_2',6);                                 % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','secant');                                                             % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.5;
options.opt_cont.step_control = 'angle';
% options.opt_cont.step_control_param = [2, 5]; 


% Continuation and postprocessing
timer = tic;                                    % Record current time
[S1,DYN1] = costar(options);                    % Calculate initial solution and continue the curve
time1 = toc(timer);                             % Display elapsed time since tic

benchmark_postprocess_quasiperiodic(DYN1,S1);   % Postprocessing


% Comparison with shooting
%{
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param);                                                   % Properties of the system
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','shooting','cont','on','stability','off',...       % Properties of the solution
                              'auto_freq',auto_freq,'non_auto_freq',non_auto_freq,'act_param',active_parameter);            % Properties of the solution
options.opt_init = costaropts('ic',[1;0]);                                                                                  % Property for initial solution (Bug in code: 'iv' needs be present, put can be empty)
options.opt_approx_method = costaropts('solver','ode45','n_char',100);                                                      % Properties of approximation method SHM
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.1,'direction',1);                                          % Properties for continuation

tic                                                                     % Record current time
[S1_SHM,DYN1_SHM] = costar(options);                                    % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot1 = costaropts('zaxis', 'max2');
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Create a new continuation plot of the continuation using FDM
opt_contplot1_SHM = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1_SHM,mu1_SHM] = S1_SHM.contplot(DYN1_SHM,opt_contplot1_SHM);         % Plot continuation 1 (using SHM) into figure of continuation plot FDM

% benchmark_postprocess_quasiperiodic(DYN1_SHM,S1_SHM)
%}


%% Continuation of excitation frequency
epsilon = 0.3;                                  % Non-linear damping coefficient
s = 0.4;                                        % Excitation amplitude

mu_limit = [1.13, 1.17];  mu0 = mu_limit(2);    % Limits of continuation diagram and mu-value epsilon0 at start of continuation

param = {epsilon, s, mu0};                      % Parameter array
active_parameter = 3;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) mu;                       % Non-autonomous frequency
Fcn = @(t,z,param) vdP_qp(t,z,param);           % Right hand side of ODE

auto_freq = 1;                                  % Initial value of autonomous frequency
C1_mat = [-0.5, 0;                              % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
          0, -1.5];
S1_mat = [0, -1.5;
          0.6, 0];


% Properties
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param);                                                           % Properties of the system
options.opt_sol  = costaropts('sol_type','qps','approx_method','fdm','cont','on','stability','off','auto_freq',auto_freq,...        % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter,'display','step-control');                 % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                               % Properties for initial solution
options.opt_approx_method = costaropts('n_int_1',50,'points_1',[-4,-3,-2,-1,0,1,2],...                                              % Properties of approximation method FDM
                                       'n_int_2',50,'points_2',[-4,-3,-2,-1,0,1,2]);                                                % Properties of approximation method FDM
% options.opt_approx_method = costaropts('n_int_1',50,'scheme_1','central','approx_order_1',6,...                                   % Properties of approximation method FDM
%                                        'n_int_2',50,'scheme_2','central','approx_order_2',6);                                     % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','secant','direction',-1);                                                  % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.25;
options.opt_cont.step_control = 'angle';
% options.opt_cont.step_control_param = [2, 5]; 


% Continuation and postprocessing
timer = tic;                                    % Record current time
[S2,DYN2] = costar(options);                    % Calculate initial solution and continue the curve
time2 = toc(timer);                             % Display elapsed time since tic

benchmark_postprocess_quasiperiodic(DYN2,S2);   % Postprocessing


% Comparison with shooting
%{
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param);                                                   % Properties of the system
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','shooting','cont','on','stability','off',...       % Properties of the solution
                              'auto_freq',auto_freq,'non_auto_freq',non_auto_freq,'act_param',active_parameter);            % Properties of the solution
options.opt_init = costaropts('ic',[1;0]);                                                                                  % Property for initial solution (Bug in code: 'iv' needs be present, put can be empty)
options.opt_approx_method = costaropts('solver','ode45','n_char',300);                                                      % Properties of approximation method SHM
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.1,'direction',-1);                                         % Properties for continuation

tic                                                                     % Record current time
[S2_SHM,DYN2_SHM] = costar(options);                                    % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot2 = costaropts('zaxis', 'max2');
[s2,mu2] = S2.contplot(DYN2,opt_contplot2);                             % Create a new continuation plot of the continuation using FDM
opt_contplot2_SHM = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s2_SHM,mu2_SHM] = S2_SHM.contplot(DYN2_SHM,opt_contplot2_SHM);          % Plot continuation 2 (using SHM) into figure of continuation plot FDM

% benchmark_postprocess_quasiperiodic(DYN2_SHM,S2_SHM);
%}