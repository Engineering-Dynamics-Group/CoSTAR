%%   Example: Parable / Fold Bifurcation (equilibrium)   %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
a = 1;                                          %
b = 1;                                          %

IC = 0.1;     mu_limit = [-5, b+1];       mu0 = mu_limit(1);

param = {mu0, a, b};
active_parameter = 1;
Fcn = @(x,param) parable(x,param);


%% Properties
options.system   = costaropts('order',0,'rhs',Fcn,'param',param,'info','Continuation of Parable Equation','dim',1);                 % Properties of the System
options.opt_sol  = costaropts('cont','on','stability','on','sol_type','equilibrium','act_param',active_parameter,'display','full'); % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                                             % Property for initial solution
options.opt_cont = costaropts('pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit);                                        % Properties for continuation

% Step control options
options.opt_cont.step_width = 0.1;
options.opt_cont.step_control = 'combination';


%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Plot
benchmark_postprocess_equilibrium(DYN,S);
