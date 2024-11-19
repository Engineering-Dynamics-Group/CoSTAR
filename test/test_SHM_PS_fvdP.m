
%% Forced Van der Pol in synchronisation showing two fold bifurcations 
%% The second fold bifurction appears on an already unstable branch and does not change the overall stability 

%% Define the parameters
mu_limit = [0,2];
epsilon = 0.5;
s = 0.9;
Omega = 1.1;

IC = [1;0];
non_auto_freq = @(mu) Omega;                                                                 %non-autonomous frequencies

sw = 0.01;
dir = -1;
mu_start = mu_limit(1,2);


%% Parameters
param = {epsilon,mu_start,Omega};               %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 2;                                                                   %Index of active parameter
Fcn = @(t,z,param)vdP_qp(t,z,param);                                        %Right-hand-side of ODE

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',2);                                                                                %Properties of the System
options.opt_sol  = costaropts('cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter,'stability','on','display','iter'); %Properties of the solution
options.opt_cont = costaropts('direction',dir,'step_width',0.01,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_control','on','plot','on'); %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);                                                                             %Properties for sol_method (e.g. Shoot)
options.opt_init = costaropts('ic',IC);
options.opt_stability = costaropts('iterate_bfp','on');

                                                     %Changes the direction of continuation
%% Continuation
timer = tic;                                    % Record current time
[S,DYN] = costar(options);                      % Calculate initial solution and continue the curve
time = toc(timer);                              % Display elapsed time since tic

benchmark_postprocess_periodic(DYN,S)
