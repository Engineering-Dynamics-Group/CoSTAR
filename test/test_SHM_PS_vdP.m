%% Example: van der Pol Oscillator (periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
IC = [2;0];                                     % Initial condition for starting solution
mu_limit = [-2,-0.1];                           % Limits of continuation diagram
auto_freq = 1;                                  % Start value for autonomous frequency

param = {mu_limit(2)};                          % Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 1;                           % Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)vdP_auto_ap(t,z,param);       % Right-hand-side of ODE


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);           % Properties of the system
options.opt_sol  = costaropts('stability','off','cont','off','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);    % Properties of the solution
options.opt_init = costaropts('ic',IC);                                             % Properties for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',10,'phase_condition','poincare');      % Properties for approximation method
options.opt_stability       = costaropts('iterate_bfp','on');                       % Properties for stability


%% Single Solution
timer = tic;                                    % Record current time
[S1,DYN1] = costar(options);                    % Calculate initial solution and continue the curve
time1 = toc(timer);                             % Display elapsed time since tic


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                               % Properties of the system
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);     % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                                 % Properties for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',25,'phase_condition','integral');                     % Properties for approximation method
options.opt_cont = costaropts('step_control','angle','direction',-1,'mu_limit',mu_limit,'step_width',0.05);             % Properties for continuation
options.opt_stability       = costaropts('iterate_bfp','on');                                                           % Properties for stability


%% Continuation
timer = tic;                                    % Record current time
[S2,DYN2] = costar(options);                    % Calculate initial solution and continue the curve
time2 = toc(timer);                             % Display elapsed time since tic

% Check the Floquet multiplicators of the user-defined Jacobian
% floquet_fsolve = S2.multipliers;              % Calculate this with Jacobian from fsolve
% [S3,DYN3] = costar(options);                  % Calculate this with user-defined Jacobian
% floquet_J = S3.multipliers;
% delta_floquet_1 = abs(floquet_fsolve(1,:) - floquet_J(1,:));
% delta_floquet_2 = abs(floquet_fsolve(2,:) - floquet_J(2,:));
% rel_floquet_1 = abs((floquet_fsolve(1,:) - floquet_J(1,:)) ./ floquet_fsolve(1,:));
% rel_floquet_2 = abs((floquet_fsolve(2,:) - floquet_J(2,:)) ./ floquet_fsolve(2,:));
%
% hold on; grid on;
% plot(S2.mu,delta_floquet_1,'Linewidth',1,'DisplayName','$|\lambda_1^{fsolve} - \lambda_1^J|$');
% plot(S2.mu,delta_floquet_2,'Linewidth',1,'LineStyle','--','DisplayName','$|\lambda_2^{fsolve} - \lambda_2^J|$');
% plot(S2.mu,rel_floquet_1,'Linewidth',1,'LineStyle','-.','DisplayName','$|(\lambda_1^{fsolve} - \lambda_1^J) / \lambda_1^{fsolve}|$');
% plot(S2.mu,rel_floquet_2,'Linewidth',1,'LineStyle','--','DisplayName','$|(\lambda_2^{fsolve} - \lambda_2^J) / \lambda_2^{fsolve}|$');
% set(gcf,'Color','w');
% ylim([1e-11 1e-5]);
% set(gca,'Fontsize',12,'YScale','log');
% xlabel('Continuation Parameter $\mu = \varepsilon$','Interpreter','Latex');
% % ylabel('$|\lambda_i^{fsolve} - \lambda_i^J|$','Interpreter','Latex');
% legend('Interpreter','Latex','FontSize',12)
% title({'Vergleich der Floquet-Multiplikatoren beim autonomen van der Pol-Schwinger',...
%        'Jacobi-Matrix von fsolve vs. eigene Jacobi-Matrix (jeweils mit zentralen Differenzen berechnet)'},'Interpreter','latex')


%% Test Postprocessing
benchmark_postprocess_periodic(DYN2,S2);
