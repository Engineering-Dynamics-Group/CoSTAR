
%% Parameters
mu_limit = [0.01,2.5];                                                                 %Limits of continuation diagram

eta = mu_limit(1,1);
Di = 0.2;
Delta = 1/3;
e = 0.25;
d3 = 0.25;
Fg = 0.3924; % 0

%% DUFFING EXAMPLE
IC = [0.1;0.1;0;0];                                                                       %Initial condition for starting solution

non_auto_freq = @(mu) mu;                                                           %Non autonomous frequency, either as function of bifurcation parameter or as a constant e.g. non_auto_freq = 2*pi

param = {eta,Di,Delta,e,d3,Fg};                                                      %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 1;                                                               %Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)laval_ps(t,z,param);                                            %Right-hand-side of ODE

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);                                                                                               %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','non_auto_freq',non_auto_freq,'sol_type','periodic','approx_method','mshm','act_param',active_parameter);       %Properties of the solution
options.opt_init = costaropts('ic',IC);
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);                                                                                                %Properties for approx_method (e.g. Shoot)
options.opt_stability       = costaropts('iterate_bfp','on');
options.opt_cont = costaropts('step_control','on','direction',1,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.1,'max_cont_step',1e4);                                                             %Properties for continuation


%% Continuation
tic
[S,DYN] = costar(options);                                                                                                                                  %Calculate initial solution and continue the curve to set limits
zeit = toc;

%% Test Postprocessing  
benchmark_postprocess_periodic(DYN,S);
