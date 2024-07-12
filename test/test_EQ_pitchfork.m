
%clear variables;

%% Parameters
gamma = 0.1;

%% Pitchfork Example
IC = 0;     mu_limit = [-1.5,1.5];       mu0 = mu_limit(1);

param = {mu0,gamma};
active_parameter = 1;
Fcn =  @(x,param)pitchfork_ap(x,param);

%% Properties
options.system              = costaropts('order',0,'rhs',Fcn,'param',param,'info','continuation of pitchfork equation','dim',1);                                                                        %Properties of the System
options.opt_sol             = costaropts('stability','on','sol_type','equilibrium','act_param',active_parameter,'cont','on');                   %Properties of the solution
options.opt_init            = costaropts('ic',IC);
options.opt_cont            = costaropts('pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit);                             %Properties for continuation
options.opt_stability       = costaropts('iterate_bfp','on');
options.opt_cont.step_width = 0.05;  
options.opt_cont.step_control = 'combination';       


%% Continuation
[S1,DYN1] = costar(options);                                                                                                    %Calculate initial solution and continue the curve to set limits

options.opt_init.ic = -1;
mu0 = mu_limit(2);
options.system.param = {mu0,gamma};
options.opt_cont.direction = -1;

[S2,DYN2] = costar(options);

% Plot both solution branches in one diagram
opt_contplot1 = costaropts('zaxis', @(z) max(z(:,1)));                  % Option structure needed to plot z (and not norm(z)) against mu
[z1,mu1] = S1.contplot(DYN1,opt_contplot1);                             % Create a new continuation plot of continuation 1
opt_contplot2 = costaropts('zaxis', @(z) max(z(:,1)), 'figure', gcf);   % Option structure needed to plot z against mu and to plot continutaion 2 into already existing figure
[z2,mu2] = S2.contplot(DYN2,opt_contplot2);                             % Plot continuation 2 into figure of continuation plot 1
xlim(mu_limit)    

% Test Postprocessing
benchmark_postprocess_equilibrium(DYN1,S1);
benchmark_postprocess_equilibrium(DYN2,S2);

% opts = costaropts('zaxis','max2');
% S2.contplot(DYN2,opts);

% figure;
% plot(S2.arclength,S2.multipliers)










