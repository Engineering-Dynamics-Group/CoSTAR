%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%                  Postprocessing                   %
%                   - contplot -                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          Welcome to the CoSTAR Examples!          %
%
% The examples provide sample code to briefly show how the toolbox can be set up and how a certain CoSTAR module can be used.
% This can be helpful if you have already used the CoSTAR toolbox. 
%
% If you are not yet familiar with CoSTAR, it is highly recommended that you work with the CoSTAR tutorials instead of the example scripts.
% The tutorials comprehensively explain certain CoSTAR modules, which is why they are the perfect starting point for CoSTAR beginners.
% There is a corresponding tutorial for each example and both of them contain the same code.
%
% This example covers the postprocessing method contplot (associated tutorial: Tutorial_Postprocessing_contplot).
% It is advised to run the sections of this example script separately and to not run the complete script.
% To do this, place the cursor in the desired section to run and click "Run Section".

% addpath(genpath('..\'))                           % Add CoSTAR to MATLAB's search path



%% Computing Solutions to work with

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
D = 0.05;     kappa = 0.3;     g = 1;               % Parameters needed for the Duffing differential equation
mu_limit = [0.01, 2.5];                             % Limits of the continuation        
eta0 = mu_limit(1);                                 % Value of continuation parameter at start of continuation
param = {kappa, D, eta0, g};                        % Parameter array
active_parameter = 3;                               % Location of continuation parameter within the array
IC = [1; 0];                                        % Initial condition (point in state space) for fsolve

% Functions
non_auto_freq = @(mu) mu;                           % Non-autonomous excitation frequency
Fcn =  @(t,z,param) duffing_ap(t,z,param);          % Right-hand side of dz/dtau = f(tau,z,kappa,D,eta,g)

% Options
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of Duffing equation');   % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','on', ...     % Properties of the solution
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);                           % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                             % Property for initial solution
options.opt_approx_method = costaropts('solver','ode45');                                                           % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit);                                                                 % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)



%% Example 1: Continuation Plot

contplot_options_1 = costaropts('zaxis','max2');                % Define the options to recreate the continuation plot
contplot_output_1  = S.contplot(DYN,contplot_options_1);        % CoSTAR postprocessing for continuation plots



%% Example 2: Maximum absolute value of the state variables

% Plot the maximum absolute value of x(theta) = z_1(theta) in green:
contplot_options_2 = costaropts('zaxis',@(z) max(abs(z(:,1))),'color','g');
contplot_output_2  = S.contplot(DYN,contplot_options_2);

% Plot the maximum absolute value of x'(theta) = z_2(theta) in cyan (using rgb code [0,0.8,0.8]) and as dashed line into the same figure:
contplot_options_3 = costaropts('zaxis',@(z) max(abs(z(:,2))),'color',[0,0.8,0.8],'linestyle','--','figure',gcf);     
contplot_output_3  = S.contplot(DYN,contplot_options_3);



%% Example 3: Section of continuation curve using different resolution

contplot_options_4 = costaropts('zaxis','max2','index',155:170,'resolution',500);     
contplot_output_4  = S.contplot(DYN,contplot_options_4);

xlim([min(contplot_output_4.mu) max(contplot_output_4.mu)])     % Adapt the abscissa limits to the plotted curve
