%% Example: van der Pol Oscillator (periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
IC = [2;0];                                     % Initial condition for starting solution
mu_limit = [0.2,2];                             % Limits of continuation diagram
auto_freq = 1;                                  % Start value for autonomous frequency

param = {mu_limit(1)};                          % Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 1;                           % Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)vdP_auto_ap(t,z,param);       % Right-hand-side of ODE


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);           % Properties of the system
options.opt_sol  = costaropts('stability','on','cont','off','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);    % Properties of the solution
options.opt_init = costaropts('ic',IC);                                             % Properties for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',1);               % Properties for approximation method
options.opt_stability       = costaropts('iterate_bfp','on');                       % Properties for stability


%% Single Solution
timer = tic;                                    % Record current time
[S1,DYN1] = costar(options);                    % Calculate initial solution and continue the curve
time1 = toc(timer);                             % Display elapsed time since tic


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                               % Properties of the system
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);     % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                                 % Properties for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',1);                                                   % Properties for approximation method
options.opt_cont = costaropts('step_control','off','pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit);       % Properties for continuation
options.opt_stability       = costaropts('iterate_bfp','on');                                                           % Properties for stability


%% Continuation
timer = tic;                                    % Record current time
[S2,DYN2] = costar(options);                    % Calculate initial solution and continue the curve
time2 = toc(timer);                             % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_periodic(DYN2,S2);
