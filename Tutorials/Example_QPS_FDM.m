%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%            Quasi-Periodic Solutions               %
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
% This example covers the computation of quasi-periodic solutions using the Finite Difference Method (associated tutorial: Tutorial_QPS_FDM).
% It is advised to run the sections of this example script separately and to not run the complete script.
% To do this, place the cursor in the desired section to run and click "Run Section".

% addpath(genpath('..\'))                           % Add CoSTAR to MATLAB's search path



%%        Example No. 1: Duffing Oscillator        %%
%             Full Non-Autonomous System            %

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
D = 0.05;            kappa = 1.5;                   % Damping and nonlinear stiffness
f1 = 0.25;           f2 = 0.25;                     % Excitation amplitudes
ratio = sqrt(2);                                    % Ratio of excitation frequencies
mu_limit = [0.8, 2];                                % Limits of the continuation        
eta0 = mu_limit(1);                                 % Value of continuation parameter at start of continuation
param = {D, kappa, f1, f2, eta0, ratio};            % Parameter array
active_parameter = 5;                               % Location of continuation parameter within the array
C1_mat = [ 0,   0,  0;                              % Fourier coefficients used to create an initial value for fsolve
          0.4, 0.4, 0];
S1_mat = [0.5, 0.5, 0;
           0,   0,  0];               

% Functions
non_auto_freq = @(mu) [mu, ratio*mu];               % Non-autonomous excitation frequencies
Fcn =  @(t,z,param) duffing_ap_qp(t,z,param);       % Right-hand side of dz/dtau = f(tau,z,D,kappa,f1,f2,eta,ratio)

% Options
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','Continuation of Quasi-Periodic Duffing Equation');    % Properties of the system
options.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','off', ...      % Properties of the solution
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                           % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                               % Property for initial solution
options.opt_approx_method = costaropts('n_int_1',25,'scheme_1','central','approx_order_1',6,...                                     % Properties of approximation method
                                       'n_int_2',25,'scheme_2','central','approx_order_2',6);                                       % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.5,'pred','secant');                                                % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
solplot_options_1 = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',25,'index',[25,200]);     % "options" structure for the solplot function
solplot_output_1  = S.solplot(DYN,solplot_options_1);                                                           % Plot hyper-time surfaces Z1(theta1,theta2,mu) at desired index


% Change of the continuation parameter eta -> kappa %

% New parameters
eta_kappa = 1.5;                                        % Excitation frequency is now fixed
mu_limit_kappa = [1, 2];                                % New limits of the continuation
kappa0 = mu_limit_kappa(1);                             % New value of continuation parameter at start of continuation
param_kappa = {D, kappa0, f1, f2, eta_kappa, ratio};    % New parameter array 
active_parameter_kappa = 2;                             % New location of continuation parameter within the array
FDM_sol = S.s(:,142);                                   % An already calculated solution is used as initial value this time

% New function
non_auto_freq_kappa = @(mu) [eta_kappa, ratio*eta_kappa];              % New non-autonomous excitation frequency

% New options
options_kappa.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_kappa,'info','Continuation of Quasi-Periodic Duffing Equation - kappa');   % Properties of the system
options_kappa.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','off', ...    % Properties of the solution
                                   'non_auto_freq',non_auto_freq_kappa,'act_param',active_parameter_kappa);                             % Properties of the solution
options_kappa.opt_init = costaropts('fdm_sol',FDM_sol,'n_int_1_fdm_sol',25,'n_int_2_fdm_sol',25);                                       % Property for initial solution
options_kappa.opt_approx_method = costaropts('n_int_1',30,'points_1',[-3,-2,-1,0,1,2,3],...                                             % Properties of approximation method
                                             'n_int_2',30,'points_2',[-3,-2,-1,0,1,2,3]);                                               % Properties of approximation method
options_kappa.opt_cont = costaropts('mu_limit',mu_limit_kappa,'step_width',0.1,'pred','secant');                                        % Properties for continuation

% Continuation
[S_kappa,DYN_kappa] = costar(options_kappa);        % CoSTAR is called by costar(options)



%%        Ex. No 2: Jeffcott(-Laval) rotor         %%
%          Mixed / Half-Autonomous System           %

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
Di = 0.2;       delta = 1/3;        d3 = 0.25;      % Damping parameters
e = 0.25;       Fg = 0.3924;                        % Eccentricity and gravitational influence
mu_limit = [1.72, 2.5];                             % Limits of continuation
eta0 = mu_limit(2);                                 % Value of continuation parameter at start of continuation
param = {eta0, Di, delta, e, d3, Fg};               % Parameter array
active_parameter = 1;                               % Location of continuation parameter within the array
auto_freq = 1;                                      % Initial value of the autonomous frequency
C0 = [-0.2; -0.2; 0; 0];                            % Fourier coefficient C_0
C1_mat = [  0,    1.35, 0;                          % Fourier coefficients used to create an initial value for fsolve
          -0.35,   0,   0;
          -0.95,   0,   0;
            0,   -1.35, 0];
S1_mat = [-0.38,   0,   0;
            0,   -1.35, 0;  
            0,   -1.35, 0;
          0.875,   0,   0];

% Functions
non_auto_freq = @(mu) mu;                           % Non-autonomous excitation frequency
Fcn =  @(t,z,param) laval_qp(t,z,param);            % Right-hand side of dz/dtau = f(tau,z,eta,Di,delta,e,d3,Fg)

% Options
options.system = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param,'info','Continuation of Jeffcott-Laval Rotor');               % Properties of the system
options.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','off', ...      % Properties of the solution
                             'non_auto_freq',non_auto_freq,'auto_freq',auto_freq,'act_param',active_parameter);                     % Properties of the solution
options.opt_init = costaropts('c0',C0,'c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                       % Properties for initial solution
options.opt_approx_method = costaropts('n_int_1',30,'scheme_1','central','approx_order_1',6,...                                     % Properties of approximation method
                                       'n_int_2',30,'scheme_2','central','approx_order_2',6);                                       % Properties of the approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.5,'pred','secant','direction',-1);                                 % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
solplot_options_2 = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',30,'index',[5,22]);   % "options" structure for the solplot function
solplot_output_2  = S.solplot(DYN,solplot_options_2);                                                       % Plot hyper-time surfaces Z1(theta1,theta2,mu) at desired index


%%        Ex. No 3: van der Pol Oscillator         %%
%              Full Autonomous System               %

clear variables; clc; close all;                    % clear workspace; clear command window; close all figures

% Parameters
alpha = 0.1;             beta = 1.1;                % Parameters of the system
mu_limit = [0.1, 1.25];                             % Limits of continuation
epsilon0 = mu_limit(1);                             % Value of continuation parameter at start of continuation       
param = {epsilon0, alpha, beta};                    % Parameter array
active_parameter = 1;                               % Location of continuation parameter within the array
auto_freq = [1.04, 1.5];                            % Initial value of the autonomous frequencies
C1_mat = [0,   0,   0;                              % Fourier coefficients used to create an initial value for fsolve       
          0, -1.4,  0;
          2,   0,   0;
          0,  2.1,  0];
S1_mat = [2,   0,   0;
          0,  1.4,  0;
          0,   0,   0;
          0,  2.1,  0];            

% Function
Fcn = @(t,z,param) coupledvdp(t,z,param);           % Right-hand side of dz/dtau = f(z,epsilon,alpha,beta)

% Options
options.system   = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param,'info','Continuation of Coupled van der Pol Oscillator');   % Properties of the system
options.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','finite-difference','cont','on','stability','off', ...      % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter);                                                   % Properties of the solution
options.opt_init = costaropts('c1_matrix',C1_mat,'s1_matrix',S1_mat);                                                               % Property for initial solution
options.opt_approx_method = costaropts('n_int_1',35,'points_1',[-4,-3,-2,-1,0,1,2],...                                              % Properties of approximation method
                                       'n_int_2',35,'points_2',[-4,-3,-2,-1,0,1,2]);                                                % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',1.5,'step_control','off','pred','secant');                           % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
solplot_options_3 = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',35,'index',[5,20]);   % "options" structure for the solplot function
solplot_output_3  = S.solplot(DYN,solplot_options_3);                                                       % Plot hyper-time surfaces Z1(theta1,theta2,mu) at desired index
