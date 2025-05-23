%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%                  Postprocessing                   %
%                    - solget -                     %
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
% This example covers the postprocessing method solget (associated tutorial: Tutorial_Postprocessing_solget).
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
options.opt_approx_method = costaropts('solver','ode45','n_shoot',2);                                               % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.05);                                               % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)



%% Example 1: Solutions of all state variables

solget_options_1 = costaropts('space','hypertime','eval','all','index',100);        % Define the options
solget_output_1  = S.solget(DYN,solget_options_1);                                  % CoSTAR postprocessing for individual solutions



%% Example 2: Solutions parametrised in time (including trajectories)

% Return the first state variable z_1 = x of two different solutions (at index 30 and at index 100) with respect to time t for t in [0,25]
solget_options_2 = costaropts('space','time','eval',@(z) z(:,1),'index',[30,100],'interval',[0,25],'resolution',500);
solget_output_2  = S.solget(DYN,solget_options_2);


% Obtain the trajectory data of the solutions at index 30 and at index 100 for their respective periods:
solget_options_3 = costaropts('space','time','eval',@(z) z(:,1:2),'index',[30,100],'resolution',500);
solget_output_3  = S.solget(DYN,solget_options_3);



%% Example 3: Frequency content of a solution

% Compute the end time of the intervals
res = 2^13;                                                                         % resolution
T_30 = 2*pi/S.freq(30);                 T_100 = 2*pi/S.freq(100);                   % periods of the solutions
Delta_t_30 = 100*T_30 / res;            Delta_t_100 = 100*T_100 / res;              % time steps between two consecutive points
int_end_30 = 100*T_30 - Delta_t_30;     int_end_100 = 100*T_100 - Delta_t_100;      % end time of the intervals

% Compute the frequency content of the solution x(t) at index 30
solget_options_4 = costaropts('space','frequency','eval',@(z) z(:,1),'index',30,'interval',[0 int_end_30],'resolution',res);
solget_output_4  = S.solget(DYN,solget_options_4);

% Compute the frequency content of the solution x(t) at index 100
solget_options_5 = costaropts('space','frequency','eval',@(z) z(:,1),'index',100,'interval',[0 int_end_100],'resolution',res);
solget_output_5  = S.solget(DYN,solget_options_5);
