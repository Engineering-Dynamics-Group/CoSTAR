clear all; close all; clc;
addpath(genpath(cd));

%% Parameters
kappa = 0.3;                                    % Coefficient of non-linear stiffness
D = 0.05;                                       % Damping factor
g = 1;                                          % Amplitude of excitation

mu_limit = [0.01, 2.5];   eta0 = mu_limit(1);   % Limits of continuation diagram and mu-value at start of continuation

param = {kappa, D, eta0, g};                    % Parameter array
active_parameter = 3;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) mu;                       % Non autonomous frequency
Fcn = @(t,z,param) duffing_ap(t,z,param);       % Right-hand side of ODE

C1 = [g;0];     S1 = [0;-eta0*g];               % Fourier-coefficients to create an initial value for fsolve to find the first point on the curve


%% Properties
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','continuation of Duffing equation');           % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','on', ...    % Properties of the solution
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                   % Properties of the solution
options.opt_init = costaropts('c1',C1,'s1',S1);                                                                             % Property for initial solution
options.opt_approx_method = costaropts('n_int',200,'scheme','central','approx_order',6);                                    % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','secant','display','step_control_info');                           % Properties for continuation
options.opt_stability = costaropts('iterate_bfp','off','solver','ode45');                                                   % Properties for stability

% Step control options
options.opt_cont.step_width = 0.2;
options.opt_cont.step_width_limit = options.opt_cont.step_width .* [0.2 10];
options.opt_cont.step_control_param = [3, 7.5/180*pi]; 


%% Continuation
tic                                             % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
toc                                             % Display elapsed time since tic








