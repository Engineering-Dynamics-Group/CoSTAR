function benchmark_postprocess_equilibrium(DYN,S)

n = size(S.mu,2);                                                                                   % Get number of curve points
index = [round(n/4,0),round(n/2,0),round(3*n/4,0)];                                                 % Find three indecies to plot different solutions


%% Solget
opt_solget = costaropts('eval',@(z) z,'space','time','resolution',1e3,'index',1);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);

opt_solget = costaropts('eval',@(z) z(1),'space','hypertime','index',index);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);

%{
opt_solget = costaropts('eval','all','space','frequency','index',index,'resolution',228);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);
%}


%% Contplot
opt_contplot = costaropts('zaxis','max2','resolution',200,'Color','r');                             % Plot bifurcation diagram using max(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis','mean2','resolution',200,'figure',gcf,'Color','b');               % Plot bifurcation diagram using mean(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis','min2','resolution',200,'figure',gcf,'Color','g');                % Plot bifurcation diagram using min(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis',@(z) z(1),'resolution',200,'Color',[0,0,0]);                      % Plot bifurcation diagram with user-set function
[z_val,mu_val] = S.contplot(DYN,opt_contplot);


end