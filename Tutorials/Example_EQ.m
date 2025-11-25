%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%              Equilibrium Solutions                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          Welcome to the CoSTAR Examples!          %
%
% The examples provide sample code to briefly show how the toolbox can be set up and how a certain CoSTAR module can be used.
% This can be helpful if you have already used the CoSTAR toolbox. 
%
% If you are not yet familiar with CoSTAR, it is highly recommended that you work with the CoSTAR tutorials instead of the example scripts.
% The tutorials comprehensively explain certain CoSTAR modules, which is why they are the perfect starting point for CoSTAR beginners.
% There is a corresponding tutorial for each example and the code of both of them is identical.
%
% This example covers the computation of equilibrium solutions (associated tutorial: Tutorial_EQ).
% It is advised to run the sections of this example script separately and to not run the complete script.
% To do this, place the cursor in the desired section to run and click "Run Section".

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
options.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','Continuation of Parable Equation');     % Properties of the system
options.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','off','act_param',active_parameter);     % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                                 % Property for initial solution
options.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off');                                                % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
contplot_options = costaropts('zaxis', @(z) max(z(:,1)));       % Option structure needed to plot z (and not norm(z)) against mu
contplot_output  = S.contplot(DYN,contplot_options);            % Create a new continuation plot



%%          Ex. No 2: Pitchfork Bifurcation        %%

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
gamma = 0;                                          % Parameters needed for the pitchfork function (apart from continuation parameter mu)
IC = 0;                                             % Initial point for fsolve to find the first point on the curve
mu_limit = [-1.5, 1.5];     mu0 = mu_limit(1);      % Limits and value of continuation parameter at begin of the continuation
param = {mu0, gamma};       active_parameter = 1;   % Parameter array and location of continuation parameter within the array

% Function
Fcn =  @(z,param) pitchfork(z,param);               % Right hand side of 0 = f(z,mu)

% Options for continuation 1
options1.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','Continuation of Pitchfork Bifurcation - Part 1');  % Properties of the system
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
options2.system   = costaropts('order',0,'dim',1,'rhs',Fcn,'param',param,'info','Continuation of Pitchfork Bifurcation - Part 2');  
options2.opt_sol  = costaropts('sol_type','equilibrium','cont','on','stability','on','act_param',active_parameter);
options2.opt_init = costaropts('ic',IC);
options2.opt_cont = costaropts('mu_limit',mu_limit,'step_control','off','step_width',0.05,'direction',-1);

% Continuation 2
[S2,DYN2] = costar(options2);                       % Calling CoSTAR and performing the second continuation

% Postprocessing
contplot_options_1 = costaropts('zaxis', @(z) max(z(:,1)));                 % Option structure needed to plot z (and not norm(z)) against mu
contplot_output_1  = S1.contplot(DYN1,contplot_options_1);                  % Create a new continuation plot of continuation 1
contplot_options_2 = costaropts('zaxis', @(z) max(z(:,1)), 'figure', gcf);  % Option structure needed to plot z against mu and to plot continuation 2 into already existing figure
contplot_output_2  = S2.contplot(DYN2,contplot_options_2);                  % Plot continuation 2 into figure of continuation plot 1
xlim(mu_limit)                                                              % Rescale horizontal axis to limits of continuation parameter
