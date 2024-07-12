%% Example: Harmonic Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
D = 0.2;            ratio = 1/sqrt(2);          % Damping factor and ratio of excitation frequencies (Omega_2/Omega_1)
f1 = 1;             f2 = 1;                     % Amplitudes of excitations

mu_limit = [0.1, 3];    mu0 = mu_limit(1);      % Limits of continuation diagram and mu-value at start of continuation

param = {D, f1, f2, mu0, ratio};                % Parameter array
active_parameter = 4;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) [mu, ratio*mu];           % Non autonomous frequencies
Fcn = @(t,z,param) ho_ap_qp(t,z,param);         % Right hand side of ODE

C1_mat = [1, 1;                                 % Fourier coefficients to create an initial value for fsolve to find the first point on the curve
          0, 0];
S1_mat = [  0,    0;
          -0.1, -0.1];


%% Properties
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param);                                                           % Properties of the system
options.opt_sol  = costaropts('sol_type','qps','approx_method','fdm','cont','on','stability','off',...                              % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                          % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                               % Properties for initial solution
options.opt_approx_method = costaropts('n_int_1',20,'scheme_1','central','approx_order_1',6,...                                     % Properties of approximation method FDM
                                       'n_int_2',20,'scheme_2','central','approx_order_2',6);                                       % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','cubic','display','step_control_info');                                    % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.5;
% options.opt_cont.step_control = 'angle';
options.opt_cont.step_control_param = [3, 4/180*pi];


%% Continuation
tic                                                                     % Record current time
[S,DYN] = costar(options);                                              % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic


%% Comparison with FGM
%{
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                   % Properties of the System
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','fgm','cont','on','stability','off',...       % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                  % Properties of the solution
options.opt_init = costaropts('hmatrix',[0,1,0;0,0,1],'c0',zeros(2,1),'cmatrix',C1_mat,'smatrix',S1_mat);         % Property for initial solution
options.opt_approx_method = costaropts('n_FFT',2^6,'error_control','on','error_limit',[1e-4,1e-3],'ec_iter_max',100,'n_hh_max',100);  
options.opt_cont = costaropts('mu_limit',mu_limit);                     % Properties for continuation  ,'step_width',0.25,'step_control_param',[3, 4/180*pi]

tic                                                                     % Record current time
[S1,DYN1] = costar(options);                                            % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Plot continuation 1 (using SHM) into figure of continuation plot FDM

benchmark_postprocess_quasiperiodic(DYN1,S1)
%}


%% Comparison with shooting
%{
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                   % Properties of the System
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','shooting','cont','on','stability','off',...       % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                  % Properties of the solution
options.opt_init = costaropts('ic',[1;1]);                                                                                  % Property for initial solution (Bug in code: 'iv' needs be present, put can be empty)
options.opt_approx_method = costaropts('solver','ode45','n_char',25);                                                       % Properties for approx_method (e.g. FDM)
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.25,'step_control_param',[3, 4/180*pi]);                     % Properties for continuation

tic                                                                     % Record current time
[S2,DYN2] = costar(options);                                            % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s2,mu2] = S2.contplot(DYN2,opt_contplot1);                             % Plot continuation 1 (using SHM) into figure of continuation plot FDM

benchmark_postprocess_quasiperiodic(DYN2,S2)
%}


%% Plots and Postprocessing

% Bifurcation diagram with max(z(:,1)) and max(z(:,2))
%{
opts = costaropts('zaxis','max2');
S.contplot(DYN,opts);
opts = costaropts('zaxis',@(z)max(z(:,1)),'figure',gcf,'color','g');
S.contplot(DYN,opts);
opts = costaropts('zaxis',@(z)max(z(:,2)),'figure',gcf,'color','r');
S.contplot(DYN,opts);
%}

% Comparison between elements of S.s and S.solget
%{
idx = 30;
sol_z1 = [S.s(1:2:end,idx); S.s(1,idx)];                                % Ortskoordinaten
sol_z2 = [S.s(2:2:end,idx); S.s(2,idx)];                                % Geschwindigkeitskoordinaten
theta = linspace(0, 2*pi, DYN.opt_approx_method.n_int+1);
opts1 = costaropts('zaxis',@(z)z(:,1),'space','hypertime','index',idx);
[~,~,~,~] = S.solplot(DYN,opts1);
hold on
plot(theta, sol_z1, 'r--')
opts2 = costaropts('zaxis',@(z)z(:,2),'space','hypertime','index',idx);
[~,~,~,~] = S.solplot(DYN,opts2);
plot(theta, sol_z2, 'r--')
%}

benchmark_postprocess_quasiperiodic(DYN,S)
