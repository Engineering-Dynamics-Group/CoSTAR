%%    Example: Duffing Oscillator (periodic)   %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


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
% load('workspace_test_duffing_PS_FDM');    options.opt_init = costaropts('fdm_sol',s0);            % Use already calculated FDM solution vector as initial value
options.opt_approx_method = costaropts('n_int',200,'scheme','central','approx_order',6);                                    % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'pred','secant','display','step_control_info');                           % Properties for continuation
options.opt_stability = costaropts('iterate_bfp','off','solver','ode45');                                                   % Properties for stability

% Step control options
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.2;
options.opt_cont.step_width_limit = options.opt_cont.step_width .* [0.2 10];
% options.opt_cont.step_control = 'angle';
options.opt_cont.step_control_param = [2, 7.5];


%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic


%% Single solution at eta = 2: Using calculated solution of continuation as initial value, but different discretisation
options.system.param = {kappa, D, 2, g};
options.opt_sol.cont = 'off';
options.opt_init = costaropts('fdm_sol',S.s(:,146));
options.opt_approx_method.n_int = 250;
[S_Single,DYN_Single] = costar(options);


%% Comparison with shooting
%{
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                           % Properties of the System
options.opt_sol  = costaropts('cont','on','non_auto_freq',non_auto_freq,'stability','on','sol_type','periodic','approx_method','shooting','act_param',active_parameter);    % Properties of the solution
options.opt_init = costaropts('ic',C1);                                                             % Property for initial solution
options.opt_approx_method = costaropts('solver','ode45');                                           % Properties for approx_method (e.g. FDM)
options.opt_cont = costaropts('mu_limit',mu_limit);                                                 % Properties for continuation

tic                                                                     % Record current time
[S1,DYN1] = costar(options);                                            % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic

opt_contplot = costaropts('zaxis', 'max2');
[s,mu] = S.contplot(DYN,opt_contplot);                                  % Create a new continuation plot of the continuation using FDM
opt_contplot1 = costaropts('zaxis', 'max2', 'color', 'r', 'linestyle', '--', 'figure', gcf);
[s1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Plot continuation 1 (using SHM) into figure of continuation plot FDM

benchmark_postprocess_periodic(DYN1,S1)
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

% Postprocessing
benchmark_postprocess_periodic(DYN,S);
