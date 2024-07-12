%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%               Periodic Solutions                  %
%               - Shooting Method -                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          Welcome to the CoSTAR Examples!          %
%
% The purpose of the examples is to give a short example on how a certain type of ...
% solution (using a particular approximation method) can be calculated in CoSTAR.
%
% The following code is identical to the code of the Periodic Solutions - Shooting Method - Tutorial.
% However, most of the explanations have been omitted in order to keep this as short as possible and to reduce it to the essentials.
% It is advised to run the sections of this script separately and to not run the complete script. ...
% To do this, place the cursor in the desired section to run and click "Run Section".
%
% addpath(genpath('..\'))                           % Add CoSTAR to MATLAB's search path



%%        Example No. 1: Duffing Oscillator        %%
%               Non-Autonomous System               %

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
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','continuation of Duffing equation');   % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','on', ...     % Properties of the solution
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);                           % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                             % Property for initial solution
options.opt_approx_method = costaropts('solver','ode45');                                                           % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit);                                                                 % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
opt_solplot = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[0.5,1,2]);    % "options" structure for the solplot function
[t_Duff,z1_Duff,mu_Duff] = S.solplot(DYN,opt_solplot);                          % Plot time-dependent approximate solution z(tau,mu) for desired mu-values


% Change of the continuation parameter eta -> kappa %

% New parameters
eta_kappa = 1.5;                                    % Excitation frequency is now fixed
mu_limit_kappa = [0, 1];                            % New limits of the continuation
kappa0 = mu_limit_kappa(1);                         % New value of continuation parameter at start of continuation
param_kappa = {kappa0, D, eta_kappa, g};            % New parameter array
active_parameter_kappa = 1;                         % New location of continuation parameter within the array
C1_kappa = [-g; 0];  S1_kappa = [0; g*eta_kappa];   % New Fourier coefficients used to create an initial value for fsolve

% New function
non_auto_freq_kappa = @(mu) eta_kappa;              % New non-autonomous excitation frequency

% New options
options_kappa.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_kappa,'info','continuation of Duffing equation - kappa');   % Properties of the system
options_kappa.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','on', ...                   % Properties of the solution
                                   'non_auto_freq',non_auto_freq_kappa,'act_param',active_parameter_kappa);                             % Properties of the solution
options_kappa.opt_init = costaropts('ic',IC);                                                                                           % Properties for initial solution
options_kappa.opt_approx_method = costaropts();                                                                                         % Properties of approximation method
options_kappa.opt_cont = costaropts('mu_limit',mu_limit_kappa);                                                                         % Properties for continuation

% Continuation
[S_kappa,DYN_kappa] = costar(options_kappa);        % CoSTAR is called by costar(options)



%%        Ex. No 2: van der Pol Oscillator         %%
%                 Autonomous System                 %

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
mu_limit = [0.1, 2.5];   epsilon0 = mu_limit(1);    % Limits and value of continuation parameter at start of the continuation        
param = {epsilon0};      active_parameter = 1;      % Parameter array and location of continuation parameter within the array
auto_freq = 1;                                      % Initial value of the autonomous frequency
IC = [2; 0];                                        % Initial condition (point in state space) for fsolve

% Function
Fcn = @(t,z,param) vdP_auto_ap(t,z,param);          % Right-hand side of dz/dtau = f(z,epsilon)

% Options
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','continuation of van der Pol oscillator');   % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','shooting','cont','on','stability','off', ...            % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter);                                           % Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                                     % Property for initial solution
options.opt_approx_method = costaropts('solver','ode45');                                                                   % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit);                                                                         % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
opt_solplot = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[1,2]);    % "options" structure for the solplot function
[t_vdP,z1_vdP,mu_vdP] = S.solplot(DYN,opt_solplot);                         % Plot time-dependent approximate solution z(tau,mu) for desired mu-values
