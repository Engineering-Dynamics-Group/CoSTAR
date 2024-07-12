function benchmark_postprocess_periodic(DYN,S)

n = size(S.mu,2);                                                                                   % Get number of curve points
index = [round(n/4,0),round(n/2,0),round(3*n/4,0)];                                                 % Get three indecies to plot different solutions


%% Solget
opt_solget = costaropts('eval',@(z)z(:,1),'space','time','resolution',1e3,'index',1);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);

opt_solget = costaropts('eval',@(z) z(:,1),'space','hypertime','index',index);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);


%% Contplot
opt_contplot = costaropts('zaxis','max2','resolution',200,'Color','r');                             % Plot bifurcation diagram using max(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis','mean2','resolution',200,'figure',gcf,'Color','b');               % Plot bifurcation diagram using mean(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis','min2','resolution',200,'figure',gcf,'Color','g');                % Plot bifurcation diagram using min(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis',@(z)max(z(:,2)),'resolution',200,'Color',[0,0,0]);                % Plot bifurcation diagram with user-set function
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

%{
opt_contplot = costaropts('zaxis','max2','resolution',200);                % Plot bifurcation diagram with user-set function
[z_val,mu_val] = S.contplot(DYN,opt_contplot);
%}


%% Solplot
opt_solplot = costaropts('zaxis',@(z) z(:,1),'space','hypertime','resolution',200,'index',index(3),'Color','r');        % Plot solution z_1 in hyper-time
[hypertime_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

%{
opt_solplot = costaropts('zaxis','all','space','hypertime','resolution',200,'index',index(2:3));
[hypertime_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);
%}


opt_solplot = costaropts('zaxis',@(z) z(:,1),'space','time','resolution',100,'index',index);                            % Plot time solution z_1 to specific index
[t_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,1),'space','time','index',index(2));                                          % Plot interval of time solution z_1 to specific index
[t_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);


opt_solplot = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','index',index);                                                   % Plot trajectory to specific index
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','mu',S.mu(1,index(2)),'figure',gcf);   % Plot trajectory to specific mu-value
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z)z(:,2),'xaxis',@(z)z(:,1),'space','trajectory','resolution',104,'index',index(2),'Color','r');                     % Plot trajectory to specific index and resolution
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','index',index(2),'interval',[0,2]);                               % Plot interval of trajectory to specific index
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);


opt_solplot = costaropts('zaxis','all','space','frequency','resolution',102,'index',index,'Color','r');     % Plot frequency spectrum to specific index
[f_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis','all','space','frequency','resolution',211,'index',index);                 % Plot frequency spectrum to specific index
[f_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);


end