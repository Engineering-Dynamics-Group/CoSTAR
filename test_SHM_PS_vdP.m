clear all; close all; clc;
addpath(genpath(cd));

%% Van-der-Pol EXAMPLE
IC = [1;0];                                                                             %Initial condition for starting solution
mu_limit = [0.2,2];                                                                     %Limits of continuation diagram
auto_freq = 1;                                                                          %Start value for autonomous frequency

param = {mu_limit(2)};                                                                    %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
dir = -1;
active_parameter = 1;                                                                   %Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)vdP_auto_ap(t,z,param);                                               %Right-hand-side of ODE

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                               %Properties of the System
options.opt_sol  = costaropts('stability','off','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);               %Properties of the solution
options.opt_cont = costaropts('step_control','off','pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'direction',dir);                                                                    %Properties for continuation
options.opt_init = costaropts('ic',IC);                                                                 %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_shoot',3);                                                                                                                 %Properties for approx_method (e.g. Shoot)
options.opt_stability       = costaropts('iterate_bfp','off');                                                                                                                      %Changes the direction of continuation (uncomment only if algorithm doesn't start properly)

%% Continuation
tic
[S,DYN] = costar(options);                                                                                                                                 %Calculate initial solution and continue the curve to set limits
zeit = toc;
