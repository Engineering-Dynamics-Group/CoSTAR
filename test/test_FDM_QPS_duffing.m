%% Example: Duffing Oscillator (quasi-periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
D = 0.05;            ratio = sqrt(2);           % Damping factor and ratio of excitation frequencies (Omega_2/Omega_1)
f1 = 0.25;           f2 = 0.25;                 % Amplitudes of excitation
kappa = 1.5;                                    % Coefficient of nonlinear stiffness

mu_limit = [0.8, 2];   mu0 = mu_limit(2);       % Limits of continuation diagram and mu-value at start of continuation

% C1_mat = [ 0,   0,  0;                        % Fourier coefficients to create an initial value for fsolve to find the first point on the curve at mu = 0.8
%           0.4, 0.4, 0];
% S1_mat = [0.5, 0.5, 0;
%            0,   0,  0];
C1_mat = [  0,     0,   0;                      % Fourier coefficients to create an initial value for fsolve to find the first point on the curve at mu = 2
          -0.2, -0.07, 0];
S1_mat = [-0.1, -0.025, 0;
            0,     0,   0];

% Difficult Duffing Oscillator (same parameters as "FGM_duffing_QPS_hard and "test_Shooting_duffing_QPS")
%{
D = 0.1;            ratio = 1/sqrt(2);          % Damping factor and ratio of excitation frequencies (Omega_2/Omega_1)
f1 = 5;             f2 = 5;                     % Amplitudes of excitations
kappa = 0.2;                                    % Coefficient of nonlinear stiffness

mu_limit = [2.2, 5];   mu0 = mu_limit(2);       % Limits of continuation diagram and mu-value at start of continuation  ->  'direction',-1
                                                % In order to decrease mu_limit(1): Discretization must be refined
C1_mat = [ 0,   0,  0;                          % Fourier coefficients to create an initial value for fsolve to find the first point on the curve at mu = 5
          -1, -1.5, 0];
S1_mat = [-0.2, -0.4, 0;
            0,    0,  0];
%}

param = {D, kappa, f1, f2, mu0, ratio};         % Parameter array
active_parameter = 5;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) [mu, ratio*mu];           % Non autonomous frequencies
Fcn = @(t,z,param) duffing_ap_qp(t,z,param);    % Right-hand side of ODE


%% Properties
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of Quasi-Periodic Duffing');           % Properties of the system
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','off',...      % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                          % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                               % Properties for initial solution
% load('workspace_test_duffing_QPS_FDM');    options.opt_init = costaropts('fdm_sol',s0);           % Use already calculated FDM solution vector as initial value
options.opt_approx_method = costaropts('n_int_1',25,'scheme_1','central','approx_order_1',6,...                                     % Properties of approximation method FDM
                                       'n_int_2',25,'scheme_2','central','approx_order_2',6);                                       % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','cubic','direction',-1,'display','step_control_info');                     % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.5;                  % Difficult Duffing: step_width = 1
options.opt_cont.step_control = 'angle';
% options.opt_cont.step_control_param = [2, 3];     % Difficult Duffing: step_control_param = [2, 5]
 

%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Single solution at eta = 1.1: Using calculated solution of continuation as initial value, but different discretisation
options.system.param = {D, kappa, f1, f2, 1.1, ratio};
options.opt_sol.cont = 'off';
options.opt_init = costaropts('fdm_sol',S.s(:,129),'n_int_1_fdm_sol',25,'n_int_2_fdm_sol',25);
options.opt_approx_method.n_int_1 = 50;     options.opt_approx_method.n_int_2 = 50;
[S_Single,DYN_Single] = costar(options);


%% Comparison with FGM
%{
c_max = [0.1188,   0.1443;  -0.0104, -0.0075];      c_max2 = [c_max, 0.1.*ones(2,14)];
s_max = [-0.1049, -0.0539;  -0.0119, -0.0205];      s_max2 = [s_max, 0.1.*ones(2,14)];

K3 = [0 1 0 3 0 -1 1 -2 2 -3  -2   2 -4  -5  -4  -5  3  3  4  5  7;...
      0 0 1 0 3  2 2  1 1  2   3   3  3   2   1   6  2  4  1  0  0];

options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                           % Properties of the System
options.opt_sol  = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off',...       % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                          % Properties of the solution
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',c_max2,'Smatrix',s_max2,'Hmatrix',K3);                                      % Properties for the initial solution
options.opt_approx_method = costaropts('n_FFT',2^5,'error_control',1,'error_limit',[0,1e-2]);                                       % Properties for approx_method | Difficult Duffing: error_limit = [0,0.5]
options.opt_cont = costaropts('mu_limit',mu_limit,'maxcontstep',1e4,'direction',1);                                                 % Properties for continuation
options.opt_cont.step_width = 0.1;                                                                                                  % Difficult Duffing: step_width = 0.5

tic                                                                     % Record current time
[S1,DYN1] = costar(options);                                            % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Plot continuation 1 (using SHM) into figure of continuation plot FDM

benchmark_postprocess_quasiperiodic(DYN1,S1)
%}


%% Plots and Postprocessing
% figure()                                                              % Open new figure
% plot(1:size(S.step_width,2), S.step_width)                            % Plot step widths

%{
hold on
plot(S.mu,S.fsolve_it,'r')

plot(S.mu,S.fval)
xlabel('mu')
ylabel('norm(fval)^2')
annotation('textbox',[.5 .5 .3 .3],'String',{'mean(norm(fval)^2) = 1.4496e-10','points on curve: 311','computing time: 4.3s'},'FitBoxToText','on')
%}

benchmark_postprocess_quasiperiodic(DYN,S);
