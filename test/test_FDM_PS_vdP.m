%%  Example: van der Pol Oscillator (periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


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
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','off', ...   % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter);                                           % Properties of the solution
options.opt_init = costaropts('c0',zeros(2,1),'c1',C1,'s1',S1);                                                             % Property for initial solution
options.opt_approx_method = costaropts('n_int',200,'points',[-4,-3,-2,-1,0,1,2]);                                           % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'display','step_control_info');                                           % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
% options.opt_cont.step_width = 0.1;
% options.opt_cont.step_control = 'off';
% options.opt_cont.step_control_param = [2, 3]; 


%% Continuation
tic                                                                     % Record current time
[S,DYN] = costar(options);                                              % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic


%% Comparison with shooting

%{
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                           % Properties of the System
options.opt_sol  = costaropts('cont','on','auto_freq',auto_freq,'stability','off','sol_type','periodic','approx_method','shooting','act_param',active_parameter);    % Properties of the solution
options.opt_init = costaropts('ic',C1);                                                             % Property for initial solution
options.opt_approx_method = costaropts('solver','ode45');                                           % Properties for approx_method (e.g. FDM)
options.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off');                            % Properties for continuation

tic                                                                     % Record current time
[S1,DYN1] = costar(options);                                            % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Plot continuation 1 (using SHM) into figure of continuation plot FDM

benchmark_postprocess_periodic(DYN1,S1)
%}


%% Postprocessing
benchmark_postprocess_periodic(DYN,S);

