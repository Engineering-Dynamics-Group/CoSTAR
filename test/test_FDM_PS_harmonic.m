%%   Example: Harmonic Oscillator (periodic)   %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
D = 0.2;            f = 1;                      % Damping factor and amplitude of excitation

mu_limit = [0.1, 3];    mu0 = mu_limit(1);      % Limits of continuation diagram and mu-value at start of continuation

param = {D, f, mu0};                            % Parameter array
active_parameter = 3;                           % Defines where the continuation parameter is located within param

non_auto_freq = @(mu) mu;                       % Non autonomous frequency
Fcn = @(t,z,param) ho_ap(t,z,param);            % Right-hand side of ODE

C1 = [1;0];     S1 = [0;-0.1];                  % Fourier-coefficients to create an initial value for fsolve to find the first point on the curve


%% Properties
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param);                                                   % Properties of the system  
options.opt_sol  = costaropts('sol_type','ps','approx_method','fdm','cont','on','stability','off',...                       % Properties of the solution
                              'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                  % Properties of the solution 
options.opt_init = costaropts('c1',C1,'s1',S1);                                                                             % Properties for initial solution             
options.opt_approx_method = costaropts('n_int',100,'scheme','forward','approx_order',2);                                    % Properties of approximation method FDM
options.opt_cont = costaropts('mu_limit',mu_limit,'display','step_control_info');                                           % Properties for continuation

% Step control options  
% Available step control methods: 'off', 'on', 'corrector_iterations', 'norm_corrector', 'combination', 'angle', ('pid')
options.opt_cont.step_width = 0.25;
% options.opt_cont.step_control = 'angle';
% options.opt_cont.step_control_param = [3, 3/180*pi];                                                                                                                   


%% Continuation
tic                                                                     % Record current time
[S,DYN] = costar(options);                                              % Calculate initial solution and continue the curve
toc                                                                     % Display elapsed time since tic


%% Comparison with shooting

%{
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                           % Properties of the System
options.opt_sol  = costaropts('cont','on','non_auto_freq',non_auto_freq,'stability','off','sol_type','periodic','approx_method','shooting','act_param',active_parameter);    % Properties of the solution
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

% Bifurcation diagram with max(z(:,1)) and max(z(:,2))
%{
opts = costaropts('zaxis','max2');
S.contplot(DYN,opts);
opts = costaropts('zaxis',@(z)max(z(:,1)),'figure',gcf,'color','g');
S.contplot(DYN,opts);
opts = costaropts('zaxis',@(z)max(z(:,2)),'figure',gcf,'color','r');
S.contplot(DYN,opts);
%}

% Comparison between elements of S.s and S.solget
%{
idx = 30;
sol_z1 = [S.s(1:2:end,idx); S.s(1,idx)];                                % Ortskoordinaten
sol_z2 = [S.s(2:2:end,idx); S.s(2,idx)];                                % Geschwindigkeitskoordinaten
theta = linspace(0, 2*pi, DYN.opt_approx_method.n_int+1);
opts1 = costaropts('zaxis',@(z)z(:,1),'space','hypertime','index',idx);
[~,~,~,~] = S.solplot(DYN,opts1);
hold on
plot(theta, sol_z1, 'r--')
opts2 = costaropts('zaxis',@(z)z(:,2),'space','hypertime','index',idx);
[~,~,~,~] = S.solplot(DYN,opts2);
plot(theta, sol_z2, 'r--')
%}

% Postprocessing
benchmark_postprocess_periodic(DYN,S)
