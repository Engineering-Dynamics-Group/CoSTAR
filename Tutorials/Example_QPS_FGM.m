%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      CoSTAR                       %
%   Continuation of Solution Torus AppRoximations   %
%                                                   %
%                     Example:                      %
%            Quasi-Periodic Solutions               %
%           - Fourier-Galerkin Method -             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%          Welcome to the CoSTAR Examples!          %
%
% The purpose of the examples is to give a short example on how a certain type of ...
% solution (using a particular approximation method) can be calculated in CoSTAR.
%
% The following code is identical to the code of the Quasi-Periodic Solutions - Fourier-Galerkin Method - Tutorial.
% However, most of the explanations have been omitted in order to keep this as short as possible and to reduce it to the essentials.
% It is advised to run the sections of this script separately and to not run the complete script. ...
% To do this, place the cursor in the desired section to run and click "Run Section".
%
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
Hmatrix = [0, 1, 0, 3, 0;                           % Matrix of the harmonics to be used (NOTE: This matrix will only lead to a ...
           0, 0, 1, 0, 3];                          % good approximation if the error control is used!).
C0 = [0; 0];                                        % Fourier coefficient C_0
Cmatrix = [0.1,  0.1;                               % Fourier coefficients C_1 and C_2          
           0.08, 0.08];                             % The coefficients C_3 and C_4 corresponding to Hmatrix(:,4:5) will be guessed  
Smatrix = [ 0.1,   0.1;                             % Fourier coefficients S_1 and S_2 
           -0.08, -0.08];                           % The coefficients S_3 and S_4 corresponding to Hmatrix(:,4:5) will be guessed  

% Functions
non_auto_freq = @(mu) [mu, ratio*mu];               % Non-autonomous excitation frequencies
Fcn =  @(t,z,param) duffing_ap_qp(t,z,param);       % Right-hand side of dz/dtau = f(tau,z,D,kappa,f1,f2,eta,ratio)

% Options
options.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param,'info','continuation of quasi-periodic Duffing equation');    % Properties of the system
options.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off', ...       % Properties of the solution
                             'non_auto_freq',non_auto_freq,'act_param',active_parameter);                                           % Properties of the solution
options.opt_init = costaropts('hmatrix',Hmatrix,'c0',C0,'cmatrix',Cmatrix,'smatrix',Smatrix);                                       % Property for initial solution
options.opt_approx_method = costaropts('n_FFT',2^5,'error_control','on','error_limit',[1e-4,1e-2],'ec_iter_max',10,'n_hh_max',50);  % Properties of approximation method
options.opt_cont = costaropts('mu_limit',mu_limit);                                                                                 % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
opt_solplot = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',30,'index',[25,200]);   % "options" structure for the solplot function
[t_Duff,z1_Duff,mu_Duff] = S.solplot(DYN,opt_solplot);                                                  % Plot hyper-time surfaces Z1(theta1,theta2,mu) at desired index


% Change of the continuation parameter eta -> kappa %

% New parameters
eta_kappa = 1.5;                                        % Excitation frequency is now fixed
mu_limit_kappa = [1, 2];                                % New limits of the continuation
kappa0 = mu_limit_kappa(1);                             % New value of continuation parameter at start of continuation
param_kappa = {D, kappa0, f1, f2, eta_kappa, ratio};    % New parameter array 
active_parameter_kappa = 2;                             % New location of continuation parameter within the array
Fc0 = S.s{1};                                           % An already calculated solution is used as initial value this time
Hmatrix_kappa = S.hmatrix{1};                           % Get the corresponding Hmatrix of Fc0

% New function
non_auto_freq_kappa = @(mu) [eta_kappa, ratio*eta_kappa];              % New non-autonomous excitation frequency

% New options
options_kappa.system = costaropts('order',1,'dim',2,'rhs',Fcn,'param',param_kappa,'info','continuation of Duffing equation - kappa');       % Properties of the system
options_kappa.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off', ...         % Properties of the solution
                                   'non_auto_freq',non_auto_freq_kappa,'act_param',active_parameter_kappa);                                 % Properties of the solution
options_kappa.opt_init = costaropts('hmatrix',Hmatrix_kappa,'fc0',Fc0);                                                                     % Property for initial solution
options_kappa.opt_approx_method = costaropts('n_FFT',2^5,'error_control','on','error_limit',[1e-4,1e-2],'ec_iter_max',10,'n_hh_max',50);    % Properties of approximation method
options_kappa.opt_cont = costaropts('mu_limit',mu_limit_kappa);                                                                             % Properties for continuation

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
Hmatrix = [0, 1, 0, 3, 0;                           % Matrix of the harmonics to be used (NOTE: This matrix will only lead to a ...
           0, 0, 1, 0, 3];                          % good approximation if the error control is used!).
C0 = [-0.2; -0.2; 0; 0];                            % Fourier coefficient C_0
Cmatrix = [  0,    1.35;                            % Fourier coefficients C_1 and C_2  
           -0.35,   0;                              % The coefficients C_3 and C_4 corresponding to Hmatrix(:,4:5) will be guessed 
           -0.95,   0;
             0,   -1.35];
Smatrix = [-0.38,   0;                              % Fourier coefficients S_1 and S_2 
             0,   -1.35;                            % The coefficients S_3 and S_4 corresponding to Hmatrix(:,4:5) will be guessed  
             0,   -1.35;
           0.875,   0  ];

% Functions
non_auto_freq = @(mu) mu;                           % Non-autonomous excitation frequency
Fcn =  @(t,z,param) laval_qp(t,z,param);            % Right-hand side of dz/dtau = f(tau,z,eta,Di,delta,e,d3,Fg)

% Options
options.system = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param,'info','continuation of Jeffcott-Laval rotor');               % Properties of the system
options.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off', ...       % Properties of the solution
                             'non_auto_freq',non_auto_freq,'auto_freq',auto_freq,'act_param',active_parameter);                     % Properties of the solution
options.opt_init = costaropts('hmatrix',Hmatrix,'c0',C0,'cmatrix',Cmatrix,'smatrix',Smatrix);                                       % Property for initial solution
options.opt_approx_method = costaropts('n_FFT',2^6,'error_control','on','error_limit',[1e-4,1e-2],'ec_iter_max',10, ...             % Properties of the approximation method
                                       'n_hh_max',50,'phasecond','int_poincare');                                                   % Properties of the approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.01,'direction',-1);                                                % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
opt_solplot = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',30,'index',[5,50]);     % "options" structure for the solplot function
[t_JLR,z1_JLR,mu_JLR] = S.solplot(DYN,opt_solplot);                                                     % Plot hyper-time surfaces Z1(theta1,theta2,mu) at desired index



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
Hmatrix = [0, 1, 0, 3, 0;                           % Matrix of the harmonics to be used (NOTE: This matrix will only lead to a ...
           0, 0, 1, 0, 3];                          % good approximation if the error control is used!).
C0 = [0; 0; 0; 0];                                  % Fourier coefficient C_0
Cmatrix = [-0.2,  0.1;                              % Fourier coefficients C_1 and C_2  
           -0.1, -1;                                % The coefficients S_3 and S_4 corresponding to Hmatrix(:,4:5) will be guessed 
            2,   -0.15;
            0.2,  2.25];
Smatrix = [ 2,   -0.1;                              % Fourier coefficients S_1 and S_2  
            0.2,  1.5;                              % The coefficients S_3 and S_4 corresponding to Hmatrix(:,4:5) will be guessed 
            0.2, -0.15;
            0.1,  1.5];          

% Function
Fcn = @(t,z,param) coupledvdp(t,z,param);           % Right-hand side of dz/dtau = f(z,epsilon,alpha,beta)

% Options
options.system   = costaropts('order',1,'dim',4,'rhs',Fcn,'param',param,'info','continuation of coupled van der Pol oscillator');   % Properties of the system
options.opt_sol = costaropts('sol_type','quasiperiodic','approx_method','fourier-galerkin','cont','on','stability','off', ...       % Properties of the solution
                             'auto_freq',auto_freq,'act_param',active_parameter);                                                   % Properties of the solution
options.opt_init = costaropts('hmatrix',Hmatrix,'c0',C0,'cmatrix',Cmatrix,'smatrix',Smatrix);                                       % Property for initial solution
options.opt_approx_method = costaropts('n_FFT',2^6,'error_control','on','error_limit',[1e-3,1e-1],'ec_iter_max',10, ...             % Properties of the approximation method
                                       'n_hh_max',50,'phasecond','int_poincare');                                                   % Properties of the approximation method
options.opt_cont = costaropts('mu_limit',mu_limit,'step_width',0.01);                                                               % Properties for continuation

% Continuation
[S,DYN] = costar(options);                          % CoSTAR is called by costar(options)

% Postprocessing
opt_solplot = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',40,'index',[5,55]);     % "options" structure for the solplot function
[t_cvdP,z1_cvdP,mu_cvdP] = S.solplot(DYN,opt_solplot);                                                  % Plot hyper-time surfaces Z1(theta1,theta2,mu) at desired index
