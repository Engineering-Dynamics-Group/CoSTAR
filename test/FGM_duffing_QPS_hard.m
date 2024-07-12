clear variables; clc; close all;                                        % clear workspace; clear command window; close all figures
%addpath(genpath('D:\Daten\Uni Kassel\CoSTAR\Matlab\v2.1.2.10_J1'))      % Add main folder of CoSTAR and all subfolders to search path


%% Parameters
kappa = 0.2;
D = 0.1;
s1 = 5;
s2 = 5;
ratio = 1./sqrt(2);


%% Duffing Example
IC = 0.1.*ones(2,1);                                                    % Initial condition
mu_limit = [1,5];                                                       % Limits of continuation diagram
non_auto_freq = @(mu) [mu mu*ratio];                                    % Non-autonomous frequencies
 
param = {D,kappa,s1,s2,mu_limit(1,2),ratio};                            % Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 5;                                                   % Index of active parameter
Fcn = @(t,z,param)duffing_qp(t,z,param);                                % Right-hand-side of ODE


%% Init data
c_max =     [0.1188,    0.1443;
            -0.0104,   -0.0075];
s_max =     [-0.1049,   -0.0539;
             -0.0119,   -0.0205];

K3 = [0 1 0 3 0 -1 1 -2 2 -3  -2   2 -4  -5  -4  -5  3  3  4  5  7;...
      0 0 1 0 3  2 2  1 1  2   3   3  3   2   1   6  2  4  1  0  0];

% Higher harmonics will be guessed
c_max2 = [c_max, 0.1.*ones(2,14)];
s_max2 = [s_max, 0.1.*ones(2,14)];
    

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                           % Properties of the System
options.opt_sol  = costaropts('cont','on','stability','off','non_auto_freq',non_auto_freq,'sol_type','quasiperiodic','approx_method','fourier-galerkin','act_param',active_parameter);   % Properties of the solution
options.opt_init = costaropts('C0',zeros(2,1),'Cmatrix',c_max2,'Smatrix',s_max2,'Hmatrix',K3);                                                                      % Properties for the initial solution
options.opt_approx_method = costaropts('n_FFT',2^5,'error_control','on','error_limit',[0,0.5]);                                                                        % Properties for approx_method (e.g. Shoot)
options.opt_cont = costaropts('pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'max_cont_step',1e4,'direction',-1);                         % Properties for continuation

% Step control
options.opt_cont.step_width = 0.5;
options.opt_cont.step_control = 'angle';
% 'off', 'corrector_iterations', 'norm_corrector_predictor', 'combination', 'angle', 'pid'
%options.opt_cont.it_nominal = 2;


%% Continuation    
tic                                             % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve to set limits
toc                                             % Display elapsed time since tic

% error_limit oder n_hh können angepasst werden, um Verfolgung zu beeinflussen
% mu_limit kann angepasst werden, um den numerisch sehr aufwändigen Bereich für mu < 1.5 abzuschneiden

