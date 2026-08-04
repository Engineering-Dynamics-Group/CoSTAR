%%         Periodic Example: Langford          %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
IC = [0.0323;-1.1812;-0.1086];                  % Initial condition for starting solution
mu_limit = [0.03,0.1];                          % Limits of continuation diagram
auto_freq = 0.873336;                           % Start value for autonomous frequency

beta = 0.7;
omega = 3.5;
rho   = 0.25;
epsilon = mu_limit(1,1);

param = {beta,omega,rho,epsilon};               % Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value
active_parameter = 4;                           % Which parameter is the bifurcation parameter?
Fcn = @(t,z,param) langford(t,z,param);         % Right-hand-side of ODE


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',3);
options.opt_sol  = costaropts('stability','on','cont','off','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);
options.opt_init = costaropts('ic',IC);
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);
options.opt_stability = costaropts('iterate_bfp','on');


%% Continuation
timer = tic;                                    % Record current time
[S1,DYN1] = costar(options);                    % Calculate initial solution and continue the curve
time1 = toc(timer);                             % Display elapsed time since tic


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',3);
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);
options.opt_init = costaropts('ic',IC);
options.opt_cont = costaropts('step_control','on','predictor','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.05,'step_width_limit',[1e-5,0.01]);
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);
options.opt_stability = costaropts('iterate_bfp','on');


%% Continuation
timer = tic;                                    % Record current time
[S2,DYN2] = costar(options);                    % Calculate initial solution and continue the curve
time2 = toc(timer);                             % Display elapsed time since tic


%% Test Postprocessing
benchmark_postprocess_periodic(DYN2,S2);
