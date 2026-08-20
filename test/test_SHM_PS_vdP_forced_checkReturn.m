%%  Example: van der Pol Oscillator (forced)   %%

% Forced Van der Pol in synchronisation showing two fold bifurcations 
% The second fold bifurction appears on an already unstable branch and does not change the overall stability 

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path

%% Parameters
mu_limit = [0.5,2];
epsilon = 1;
s = 0.5;

IC = [-1.3513;1.0979];
non_auto_freq = @(mu) mu;                                                   % Non-autonomous frequencies

sw = 0.4;
dir = -1;
mu_start = 1;

param = {epsilon,s,mu_start};                                               % Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 3;                                                       % Index of active parameter
Fcn = @(t,z,param)vdP_forced(t,z,param);                                    % Right-hand-side of ODE

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                %Properties of the System
options.opt_sol  = costaropts('cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter,'stability','on','display','iter-detailed'); %Properties of the solution
options.opt_cont = costaropts('direction',dir,'step_width',sw,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_control','on','plot','on'); %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);                                                                             %Properties for sol_method (e.g. Shoot)
options.opt_init = costaropts('ic',IC);
options.opt_stability = costaropts('iterate_bfp','on');

options.opt_cont.max_cont_step = 300;

%% Continuation
timer = tic;                                                                % Record current time
[S,DYN] = costar(options);                                                  % Calculate initial solution and continue the curve
time = toc(timer);                                                          % Elapsed time


% If two many points are continued throw an error if returnCheck does not
% stop the continuation. For now set two iter = 271 (269 in console, needs
% fixing)
%
if(size(S.mu,2)>271)                                                        
    error('Return Check does not work.');
end