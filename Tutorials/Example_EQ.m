%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%              Equilibrium Solutions                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          Welcome to the CoSTAR Examples!          %
%
% The purpose of the examples is to give a short example on how a certain type of ...
% solution (using a particular approximation method) can be calculated in CoSTAR.
%
% The following code is identical to the code of the Equilibrium Solutions Tutorial.
% However, most of the explanations have been omitted in order to keep this as short as possible and to reduce it to the essentials.
% It is advised to run the sections of this script separately and to not run the complete script. ...
% To do this, place the cursor in the desired section to run and click "Run Section".
%
% addpath(genpath('..\'))                           % Add CoSTAR to MATLAB's search path



%%              Example No. 1: Parable             %%

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
a = 1;                  b = 1;                      % Parameters needed for the parable function (apart from continuation parameter mu)
IC = 0.1;                                           % Initial point for fsolve to find the first point on the curve
mu_limit = [-5, b];     mu0 = mu_limit(1);          % Limits and value of continuation parameter at begin of the continuation
param = {mu0, a, b};    active_parameter = 1;       % Parameter array and location of continuation parameter within the array

% Function
Fcn =  @(z,param) parable(z,param);                 % Right-hand side of 0 = f(z,mu)

% Options
options.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','continuation of parable equation');     % Properties of the system
options.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','off','act_param',active_parameter);     % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                                 % Property for initial solution
options.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off');                                                % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
opt_contplot = costaropts('zaxis', @(z) max(z(:,1)));   % Option structure needed to plot z (and not norm(z)) against mu
[z,mu] = S.contplot(DYN,opt_contplot);                  % Create a new continuation plot



%%          Ex. No 2: Pitchfork Bifurcation        %%

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
gamma = 0;                                          % Parameters needed for the pitchfork function (apart from continuation parameter mu)
IC = 0;                                             % Initial point for fsolve to find the first point on the curve
mu_limit = [-1.5, 1.5];     mu0 = mu_limit(1);      % Limits and value of continuation parameter at begin of the continuation
param = {mu0, gamma};       active_parameter = 1;   % Parameter array and location of continuation parameter within the array

% Function
Fcn =  @(z,param) pitchfork_ap(z,param);            % Right hand side of 0 = f(z,mu)

% Options for continuation 1
options1.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','continuation of pitchfork bifurcation - part 1');  % Properties of the system
options1.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','on','act_param',active_parameter);                 % Properties of the solution
options1.opt_init = costaropts('ic',IC);                                                                                            % Property for initial solution 
options1.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off');                                                           % Properties for continuation

% Continuation 1
[S1,DYN1] = costar(options1);                       % Calling CoSTAR and performing the first continuation

% Parameter changes for continuation 2
IC = 1;                                             % Initial point for fsolve to find the first point on the new branch
mu0 = mu_limit(2);                                  % Continuation starts at upper limit of mu this time
param = {mu0, gamma};                               % Parameter array needs to be updated

% Options for continuation 2
options2.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','continuation of pitchfork bifurcation - part 2');  
options2.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','on','act_param',active_parameter);
options2.opt_init = costaropts('ic',IC);
options2.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off','step_width',0.05,'direction',-1);

% Continuation 2
[S2,DYN2] = costar(options2);                       % Calling CoSTAR and performing the second continuation

% Postprocessing
opt_contplot1 = costaropts('zaxis', @(z) max(z(:,1)));                  % Option structure needed to plot z (and not norm(z)) against mu
[z1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Create a new continuation plot of continuation 1
opt_contplot2 = costaropts('zaxis', @(z) max(z(:,1)), 'figure', gcf);   % Option structure needed to plot z against mu and to plot continutaion 2 into already existing figure
[z2,mu2] = S2.contplot(DYN2,opt_contplot2);                             % Plot continuation 2 into figure of continuation plot 1
xlim(mu_limit)                                                          % Rescale horizontal axis to limits of continuation parameter