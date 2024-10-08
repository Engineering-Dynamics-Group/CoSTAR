%% Example: Duffing Oscillator (periodic) %%

% clear variables; clc; close all;              % clear workspace; clear command window; close all figures
% addpath(genpath('..\'))                       % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
kappa = 1.5;
D = 0.2;
g = 1;

IC = [0;1];                                                                       %Initial condition for starting solution
mu_limit = [0.1,2.5];                                                              %Limits of continuation diagram
non_auto_freq = @(mu) mu;                                                           %Non autonomous frequency, either as function of bifurcation parameter or as a constant e.g. non_auto_freq = 2*pi

param = {kappa,D,mu_limit(1),g};                                                      %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 3;                                                               %Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)duffing_ap(t,z,param);                                            %Right-hand-side of ODE


%% Properties for single solution
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                               %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','off','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);       %Properties of the solution
options.opt_init = costaropts('ic',IC);
options.opt_approx_method = costaropts('solver','ode45','n_shoot',6);                                                                                                %Properties for approx_method (e.g. Shoot)
options.opt_stability     = costaropts('iterate_bfp','on');


%% Single Solution
tic
[S,DYN] = costar(options);                                                                                                                                  %Calculate initial solution and continue the curve to set limits
zeit = toc;

% opts = struct('space','solution','eval','all');
% [s,mu] = S.solget(DYN,opts);


%% Properties for continuation
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                               %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);       %Properties of the solution
options.opt_init = costaropts('ic',IC);
options.opt_approx_method = costaropts('solver','ode45','n_shoot',6);                                                                                                %Properties for approx_method (e.g. Shoot)
options.opt_cont = costaropts('step_width',0.05,'step_width_limit',[0.01,0.1],'step_control','angle','pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit);   %Properties for continuation
options.opt_stability = costaropts('iterate_bfp','on');                                                                                                                


%% Continuation
tic
[S,DYN] = costar(options);                                                                                                                                  %Calculate initial solution and continue the curve to set limits
zeit = toc;

%     opts = struct('space','solution','eval','all');
%     [s,mu] = S.solget(DYN,opts);

% opts = costaropts('zaxis',@(z)max(z(:,1)),'color','r');
% [s,mu] = S.contplot(DYN,opts);
% 
% 
% figure;
% plot(mu,s,'r','Linewidth',2);
% set(gcf,'Color','w');
% grid on;
% xlim(mu_limit);
% set(gca,'Fontsize',12);
% xlabel('Continuation parameter $\mu$','Interpreter','Latex');
% ylabel('Amplitude $\Vert z \Vert$','Interpreter','Latex');


%% Test Postprocessing
benchmark_postprocess_periodic(DYN,S);
