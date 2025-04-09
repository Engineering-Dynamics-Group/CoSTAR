%% Example: Duffing Oscillator (periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
kappa = 0.3;        % kappa = 1.5
D = 0.05;           % D = 0.2
g = 1;

IC = [0;1];                                     % Initial condition for starting solution
mu_limit = [0.01,2.5];                          % Limits of continuation diagram | [0.1,2.5];
non_auto_freq = @(mu) mu;                       % Non autonomous frequency, either as function of bifurcation parameter or as a constant e.g. non_auto_freq = 2*pi

param = {kappa,D,mu_limit(1),g};                % Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 3;                           % Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)duffing_ap(t,z,param);        % Right-hand-side of ODE


%% Properties for single solution
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);           % Properties of the system
options.opt_sol  = costaropts('stability','off','cont','off','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);    % Properties of the solution
options.opt_init = costaropts('ic',IC);                                             % Properties for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',6);               % Properties for approximation method
options.opt_stability     = costaropts('iterate_bfp','on');                         % Properties for stability


%% Single Solution
timer = tic;                                    % Record current time
[S1,DYN1] = costar(options);                    % Calculate initial solution and continue the curve
time1 = toc(timer);                             % Display elapsed time since tic


%% Properties for continuation
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);           % Properties of the system
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);    % Properties of the solution
options.opt_init = costaropts('ic',IC);                                             % Properties for initial solution
options.opt_approx_method = costaropts('solver','ode45','n_shoot',6);               % Properties for approximation method
options.opt_cont = costaropts('step_width',0.1,'step_width_limit',[0.02,0.5],'step_control','angle','step_control_param',[2,3],'mu_limit',mu_limit);                    % Properties for continuation
options.opt_stability = costaropts('iterate_bfp','on');                             % Properties for stability                                                                                      


%% Continuation
timer = tic;                                    % Record current time
[S2,DYN2] = costar(options);                    % Calculate initial solution and continue the curve
time2 = toc(timer);                             % Display elapsed time since tic

% Check the Floquet multiplicators of the user-defined Jacobian
% floquet_fsolve = S2.multipliers;              % Calculate this with Jacobian from fsolve
% floquet_fsolve_raw = floquet_fsolve;
% floquet_fsolve = [floquet_fsolve_raw(:,1:241), zeros(2,1), floquet_fsolve_raw(:,242:297), zeros(2,1), floquet_fsolve_raw(:,298:end)];
% [S3,DYN3] = costar(options);                  % Calculate this with user-defined Jacobian
% floquet_J = S3.multipliers;    
% delta_floquet_1 = abs(floquet_fsolve(1,:) - floquet_J(1,:));
% delta_floquet_2 = abs(floquet_fsolve(2,:) - floquet_J(2,:));
% rel_floquet_1 = abs((floquet_fsolve(1,:) - floquet_J(1,:)) ./ floquet_fsolve(1,:));
% rel_floquet_2 = abs((floquet_fsolve(2,:) - floquet_J(2,:)) ./ floquet_fsolve(2,:));
%
% subplot(1,2,1)
% hold on; grid on;
% plot(1:numel(S3.n_unstable),delta_floquet_1,'Linewidth',1,'DisplayName','$\lambda_1$');
% plot(1:numel(S3.n_unstable),delta_floquet_2,'Linewidth',1,'LineStyle','--','DisplayName','$\lambda_2$');
% set(gcf,'Color','w');
% ylim([1e-9 1e-4]);
% set(gca,'Fontsize',12,'YScale','log');
% xlabel('Index','Interpreter','Latex');
% ylabel('$|\lambda_i^{fsolve} - \lambda_i^J|$','Interpreter','Latex');
% legend('Interpreter','Latex','FontSize',12)
% title({'Floquet-Multiplikatoren Duffing: Absoluter Fehler','Jacobi-Matrix von fsolve vs. eigene Jacobi-Matrix (zentrale Differenzen)'},'Interpreter','latex')
% subplot(1,2,2)
% hold on; grid on;
% plot(1:numel(S3.n_unstable),rel_floquet_1,'Linewidth',1,'DisplayName','$\lambda_1$');
% plot(1:numel(S3.n_unstable),rel_floquet_2,'Linewidth',1,'LineStyle','--','DisplayName','$$\lambda_2$');
% set(gcf,'Color','w');
% ylim([2e-9 2e-4]);
% set(gca,'Fontsize',12,'YScale','log');
% xlabel('Index','Interpreter','Latex');
% ylabel('$|(\lambda_i^{fsolve} - \lambda_i^J) / \lambda_i^{fsolve}|$','Interpreter','Latex');
% legend('Interpreter','Latex','FontSize',12)
% title({'Floquet-Multiplikatoren Duffing: Relativer Fehler','Jacobi-Matrix von fsolve vs. eigene Jacobi-Matrix (zentrale Differenzen)'},'Interpreter','latex')


%% Test Postprocessing
benchmark_postprocess_periodic(DYN2,S2);
