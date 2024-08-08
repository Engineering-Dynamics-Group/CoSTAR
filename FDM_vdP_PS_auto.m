clear all; close all; clc;
addpath(genpath(cd));

%% Parameters
mu_limit = [0.01, 3];                           % Limits of continuation diagram
epsilon0 = mu_limit(1);                         % mu-value at start of continuation
   
param = {epsilon0};                             % Parameter array
active_parameter = 1;                           % Defines where the continuation parameter is located within param

Fcn = @(t,z,param) vdP_auto_ap(t,z,param);

auto_freq = 1;                                  % Initial autonomous frequency
C1 = [2;0];   S1 = [0;-2];                      % Fourier-coefficients to create an initial value for fsolve to find the first point on the curve


%% Properties
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','continuation of van der Pol oscillator');   % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','on', ...   % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter);                                           % Properties of the solution
options.opt_init = costaropts('c0',zeros(2,1),'c1',C1,'s1',S1);                                                             % Property for initial solution
options.opt_approx_method = costaropts('n_int',200,'points',[-4,-3,-2,-1,0,1,2]);                                           % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'display','step_control_info');                                           % Properties for continuation

%% Continuation
tic                                                                     % Record current time
[S,DYN] = costar(options);                                              % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic



figure('Color','white')
plot(S.arclength,abs(S.multipliers),'LineWidth',2)
grid on

