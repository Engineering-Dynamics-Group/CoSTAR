
%% Van-der-Pol EXAMPLE
IC = [0.0323;-1.1812;-0.1086];                                                                             %Initial condition for starting solution
mu_limit = [0.03,0.1];                                                                  %Limits of continuation diagram
auto_freq = 0.873336;                                                                   %Start value for autonomous frequency


beta = 0.7;
omega = 3.5;
rho   = 0.25;
epsilon = mu_limit(1,1);


param = {beta,omega,rho,epsilon};                                                                  %Parameter vector, all constant parameters are set here, the bifurcation parameter gets its starting value (here the left corner of bifurcation diagram)
active_parameter = 4;                                                                   %Which parameter is the bifurcation parameter?
Fcn = @(t,z,param)langford(t,z,param);                                               %Right-hand-side of ODE

%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',3);                                                                                               %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','off','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);               %Properties of the solution
options.opt_init = costaropts('ic',IC);                                                                 %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);                                                                                                                 %Properties for approx_method (e.g. Shoot)
options.opt_stability = costaropts('iterate_bfp','on');                                                                                                                      %Changes the direction of continuation (uncomment only if algorithm doesn't start properly)

%% Continuation
tic
[S,DYN] = costar(options);                                                                                                                                 %Calculate initial solution and continue the curve to set limits
zeit = toc;


%% Properties
options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',3);                                                                                               %Properties of the System
options.opt_sol  = costaropts('stability','on','cont','on','auto_freq',auto_freq,'sol_type','periodic','approx_method','shm','act_param',active_parameter);               %Properties of the solution
options.opt_init = costaropts('ic',IC);
options.opt_cont = costaropts('step_control','on','pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.05,'step_width_limit',[1e-5,0.01]);                                                                    %Properties for continuation
options.opt_approx_method = costaropts('solver','ode45','n_shoot',5);                                                                                                                 %Properties for approx_method (e.g. Shoot)
options.opt_stability       = costaropts('iterate_bfp','on');                                                                                                                       %Changes the direction of continuation (uncomment only if algorithm doesn't start properly)




%% Continuation
tic
[S,DYN] = costar(options);                                                                                                                                 %Calculate initial solution and continue the curve to set limits
zeit = toc;

%% Test Postprocessing
benchmark_postprocess_periodic(DYN,S);







