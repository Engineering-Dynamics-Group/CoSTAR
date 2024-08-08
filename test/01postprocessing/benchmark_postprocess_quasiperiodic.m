function benchmark_postprocess_quasiperiodic(DYN,S)

n = size(S.mu,2);                                                                                   % Get number of curve points
index = [round(n/4,0),round(n/2,0),round(3*n/4,0)];                                                 % Find three indecies to plot different solutions


%% Solget
opt_solget = costaropts('eval',@(z)z(:,1),'space','time','resolution',1e3,'index',1);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);

opt_solget = costaropts('eval',@(z)z(:,:,1),'space','hypertime','index',index);
[z_val,mu_val,t_val,options] = S.solget(DYN,opt_solget);


%% Contplot
opt_contplot = costaropts('zaxis','max2','resolution',30,'Color','r');                              % Plot bifurcation diagram using max(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis','mean2','resolution',30,'figure',gcf,'Color','b');                % Plot bifurcation diagram using mean(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis','min2','resolution',30,'figure',gcf,'Color','g');                 % Plot bifurcation diagram using min(norm())
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

% opt_contplot = costaropts('zaxis',@(z)max(z(:,2)),'resolution',30,'Color',[0,0,0]);                 % Plot bifurcation diagram with user-set function
% [z_val,mu_val] = S.contplot(DYN,opt_contplot);

opt_contplot = costaropts('zaxis',@(z)max(z(:,:,2),[],'all'),'resolution',30,'Color',[0,0,0]);      % Plot bifurcation diagram with user-set function
[z_val,mu_val] = S.contplot(DYN,opt_contplot);

%{
opt_contplot = costaropts('zaxis','max2','resolution',50,'index','all');
[z_val,mu_val] = S.contplot(DYN,opt_contplot);
%}


%% Solplot
opt_solplot = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',20,'index',index(2:3));        % Plot solution z_1 in hyper-time
[hypertime_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,:,2),'space','hypertime','resolution',50,'index',index(2));          % Plot solution z_2 in hyper-time
[hypertime_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis','all','space','hypertime','resolution',20,'index',index(2:3));                % Plot all solutions in hyper-time
[hypertime_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

%{
opt_solplot = costaropts('zaxis','all','space','hypertime','resolution',[25,25],'index',index(2:3));                    % Plot for individual checks
[hypertime_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);
%}

opt_solplot = costaropts('zaxis',@(z)z(:,1),'space','time','resolution',100,'index',index);                                % Plot time solution z_1 to specific index
[t_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z)z(:,1),'space','time','resolution',1000,'index',1,'interval',[0,100],'Color','k');    % Plot interval of time solution z_1 to specific index and resolution
[t_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);


opt_solplot = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','index',index);                                                   % Plot trajectory to specific index
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','mu',S.mu(1,round(size(S.mu,2)/2,0)),'figure',gcf,'Color','m');   % Plot trajectory to specific mu-value
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z)z(:,2),'xaxis',@(z)z(:,1),'space','trajectory','resolution',104,'index',index(2),'Color','r');                    % Plot trajectory to specific index and resolution
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','resolution',1000,'index',index(2),'interval',[0,30]);            % Plot interval of trajectory to specific index and resolution
[x_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);


opt_solplot = costaropts('zaxis','all','space','frequency','resolution',2^(8)+1,'index',index,'Color','r');                       % Plot frequency spectrum to specific index
[f_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);

opt_solplot = costaropts('zaxis',@(z) z(:,1),'space','frequency','resolution',2^(10),'index',index(2),'interval',[0,100]);     % Plot frequency spectrum to specific index
[f_val,z_val,mu_val,empty] = S.solplot(DYN,opt_solplot);


end