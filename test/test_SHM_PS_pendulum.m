%%  Example: Mathematical Pendulum (periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
g = 9.81;                                       % Gravitational constant
l = 15;                                         % Pendulum length
D = 0;                                          % Damping parameter

phi_0 = 0.05;               phi_end = 3.1;      % Angles where to start and end the continuation
I_0 = g/l*(1-cos(phi_0));   I_end = g/l*(1-cos(phi_end));
mu_limit = [I_0, I_end];                        % Limits of continuation

param = {g,l,D,mu_limit(1)};                    % Parameter array
active_parameter = 4;                           % Defines where the continuation parameter is located within param

Fcn = @(t,z,param) math_pendulum(t,z,param);    % Right-hand-side of ODE
I = @(z,param) 1/2.*z(2,:).^2 + param{1}/param{2}.*(1 - cos(z(1,:)));   % First integral of Fcn

auto_freq = sqrt(g/l);                          % Start value for autonomous frequency
ic = [0.1; 0];                                  % Initial point in state space


%% Properties
options.system  = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'first_integral',I,'info','Continuation of Mathematical Pendulum');   % Properties of the system
options.opt_sol = costaropts('cont','on','stability','off','sol_type','periodic','approx_method','shooting', ...            % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter,'display','full');                          % Properties of the solution
options.opt_init = costaropts('ic',ic);                                                                                     % Property for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',10,'phase_condition','integral');                          % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_control','on','step_width',0.5,'pred','secant');                    % Properties for continuation

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
% options.opt_cont.step_width = 0.1;
% options.opt_cont.step_control = 'off';
% options.opt_cont.step_control_param = [2, 3]; 


%% Continuation via first integral
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Continuation via rod length
mu_limit_2 = [g, 20];                           % When starting at g, the first orbit is nearly a circle

phi_0_2 = 0.2;                                  % Angle where to fix the first integral
I_0_2 = g/mu_limit_2(1)*(1-cos(phi_0_2));       % Value of first integral

param_2 = {g,mu_limit_2(1),D,I_0_2};            % Parameter array
active_parameter_2 = 2;                         % Defines where the continuation parameter is located within param

options_2.system  = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_2,'first_integral',I,'info','Continuation of Mathematical Pendulum');  % Properties of the system
options_2.opt_sol = costaropts('cont','on','stability','off','sol_type','periodic','approx_method','shooting', ...          % Properties of the solution
                               'auto_freq',sqrt(g/mu_limit_2(1)),'act_param',active_parameter_2,'display','full');          % Properties of the solution
options_2.opt_init = costaropts('ic',ic);                                                                                   % Property for initial solution
options_2.opt_approx_method = costaropts('solver','ode45','n_shoot',10,'phase_condition','integral');                        % Properties of approximation method
options_2.opt_cont = costaropts('mu_limit',mu_limit_2,'step_control','off','step_width',1,'pred','secant');               % Properties for continuation

[S_2,DYN_2] = costar(options_2);                % Calculate initial solution and continue the curve


%% Postprocessing
benchmark_postprocess_periodic(DYN,S);
benchmark_postprocess_periodic(DYN_2,S_2);
