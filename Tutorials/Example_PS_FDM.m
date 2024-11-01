%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%               Periodic Solutions                  %
%          - Finite Difference Method -             %
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
% This example covers the computation of periodic solutions using the Finite Difference Method (associated tutorial: Tutorial_PS_FDM).
% It is advised to run the sections of this example script separately and to not run the complete script.
% To do this, place the cursor in the desired section to run and click "Run Section".

% addpath(genpath('..\'))                           % Add CoSTAR to MATLAB's search path



%%        Example No. 1: Duffing Oscillator        %%
%               Non-Autonomous System               %

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
D = 0.05;     kappa = 0.3;     g = 1;               % Parameters needed for the Duffing differential equation
mu_limit = [0.01, 2.5];                             % Limits of the continuation        
eta0 = mu_limit(1);                                 % Value of continuation parameter at begin of continuation
param = {kappa, D, eta0, g};                        % Parameter array
active_parameter = 3;                               % Location of continuation parameter within the array
C1 = [g; 0];     S1 = [0; -g*eta0];                 % Fourier coefficients used to create an initial value for fsolve

% Functions
non_auto_freq = @(mu) mu;                           % Non-autonomous excitation frequency
Fcn =  @(t,z,param) duffing_ap(t,z,param);          % Right-hand side of dz/dtau = f(tau,z,kappa,D,eta,g)

% Options
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of Duffing Equation');           % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','on', ...    % Properties of the solution
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                   % Properties of the solution
options.opt_init = costaropts('c1',C1,'s1',S1);                                                                             % Property for initial solution
options.opt_approx_method = costaropts('n_int',100,'scheme','central','approx_order',6);                                    % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.25,'step_control_param',[2,7.5]);                          % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
solplot_options_1 = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[0.5,1,2]);      % "options" structure for the solplot function
solplot_output_1  = S.solplot(DYN,solplot_options_1);                                   % Plot time-dependent approximate solution z(tau,mu) for desired mu-values


% Change of the continuation parameter eta -> kappa %

% New parameters
eta_kappa = 1.5;                                    % Excitation frequency is now fixed
mu_limit_kappa = [0, 1];                            % New limits of the continuation
kappa0 = mu_limit_kappa(1);                         % New value of continuation parameter at start of continuation
param_kappa = {kappa0, D, eta_kappa, g};            % New parameter array
active_parameter_kappa = 1;                         % New location of continuation parameter within the array
FDM_sol = S.s(:,1);                                 % An already calculated solution is used as initial value this time

% New function
non_auto_freq_kappa = @(mu) eta_kappa;              % New non-autonomous excitation frequency

% New options
options_kappa.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_kappa,'info','continuation of Duffing equation - kappa');   % Properties of the system
options_kappa.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','on', ...          % Properties of the solution
                                   'non_auto_freq',non_auto_freq_kappa,'act_param',active_parameter_kappa);                             % Properties of the solution
options_kappa.opt_init = costaropts('fdm_sol',FDM_sol);                                                                                 % Property for initial solution
options_kappa.opt_approx_method = costaropts('n_int',150,'points',[-4,-3,-2,-1,0,1,2]);                                                 % Properties of approximation method
options_kappa.opt_cont = costaropts('mu_limit',mu_limit_kappa,'step_width',0.25,'step_control_param',[2,7.5]);                          % Properties for continuation

% Continuation
[S_kappa,DYN_kappa] = costar(options_kappa);        % CoSTAR is called by costar(options)



%%        Ex. No 2: van der Pol Oscillator         %%
%                 Autonomous System                 %

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
mu_limit = [0.1, 2.5];   epsilon0 = mu_limit(1);    % Limits and value of continuation parameter at start of the continuation        
param = {epsilon0};      active_parameter = 1;      % Parameter array and location of continuation parameter within the array
auto_freq = 1;                                      % Initial value of the autonomous frequency
C1 = [2; 0];             S1 = [0; -2];              % Fourier coefficients used to create an initial value for fsolve

% Function
Fcn = @(t,z,param) vdP_auto_ap(t,z,param);          % Right-hand side of dz/dtau = f(z,epsilon)

% Options
options.system   = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of van der Pol Oscillator');   % Properties of the system
options.opt_sol = costaropts('sol_type','periodic','approx_method','finite-difference','cont','on','stability','off', ...   % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter);                                           % Properties of the solution
options.opt_init = costaropts('c1',C1,'s1',S1);                                                                             % Property for initial solution
options.opt_approx_method = costaropts('n_int',100,'points',[-4,-3,-2,-1,0,1,2]);                                           % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit);                                                                         % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
solplot_options_2 = costaropts('zaxis',@(z) z(:,1),'space','time','mu',[1,2]);      % "options" structure for the solplot function
solplot_output_2  = S.solplot(DYN,solplot_options_2);                               % Plot time-dependent approximate solution z(tau,mu) for desired mu-values
