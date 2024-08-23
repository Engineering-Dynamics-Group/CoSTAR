function benchmark_postprocess_quasiperiodic(DYN,S)

n = size(S.mu,2);                                                                                   % Get number of curve points
index = [round(n/4,0),round(n/2,0),round(3*n/4,0)];                                                 % Find three indecies to plot different solutions


%% Solget
solget_options = costaropts('eval',@(z)z(:,1),'space','time','resolution',1e3,'index',1);
solget_output  = S.solget(DYN,solget_options);

solget_options = costaropts('eval',@(z)z(:,:,1),'space','hypertime','index',index);
solget_output  = S.solget(DYN,solget_options);



%% Contplot
contplot_options = costaropts('zaxis','max2','Color','r');                          % Plot bifurcation diagram using max(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis','mean2','figure',gcf,'Color','b');            % Plot bifurcation diagram using mean(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis','min2','figure',gcf,'Color','g');             % Plot bifurcation diagram using min(norm())
contplot_output  = S.contplot(DYN,contplot_options);


contplot_options = costaropts('zaxis',@(z)max(z(:,:,1),[],'all'),'resolution',100,'Color',[0,0,0]);  % Plot bifurcation diagram with user-set function
contplot_output  = S.contplot(DYN,contplot_options);


%{
contplot_options = costaropts('zaxis','max2','resolution',50,'index','all');
contplot_output  = S.contplot(DYN,contplot_options);
%}



%% Solplot
solplot_options = costaropts('zaxis',@(z) z(:,:,1),'space','hypertime','resolution',20,'index',index(2:3));     % Plot solution z_1 in hyper-time
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,:,2),'space','hypertime','resolution',50,'index',index(2));       % Plot solution z_2 in hyper-time
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis','all','space','hypertime','resolution',20,'index',index(2:3));             % Plot all solutions in hyper-time
solplot_output  = S.solplot(DYN,solplot_options);

%{
solplot_options = costaropts('zaxis','all','space','hypertime','resolution',[25,25],'index',index(2:3));        % Plot for individual checks
solplot_output  = S.solplot(DYN,solplot_options);
%}

solplot_options = costaropts('zaxis',@(z)z(:,1),'space','time','resolution',100,'index',index);                                     % Plot time solution z_1
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z)z(:,1),'space','time','resolution',1000,'index',index(2),'interval',[0,100],'Color','k');  % Plot interval of time solution z_1 of specific index
solplot_output  = S.solplot(DYN,solplot_options);


solplot_options = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','index',index);                                                                   % Plot trajectory of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','mu',S.mu(1,index(2)),'figure',gcf);                                              % Plot trajectory of specific mu-value
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z)z(:,2),'xaxis',@(z)z(:,1),'space','trajectory','resolution',1000,'index',index(2),'interval',[0,30]);                              % Plot interval of trajectory of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','resolution',104,'index',index(2),'Color','r','linestyle','--','figure',gcf);     % Plot trajectory of specific index
solplot_output  = S.solplot(DYN,solplot_options);


solplot_options = costaropts('zaxis','all','space','frequency','index',index(2),'resolution',2^10,'interval',[0,10*2*pi/S.freq(min(index(2)))]);                    % Plot frequency spectrum of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,1),'space','frequency','index',index(1),'resolution',2^10,'interval',[0,10*2*pi/S.freq(min(index(1)))]);              % Plot frequency spectrum of z_1 of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,1),'space','frequency','index',index(3),'resolution',2^10,'interval',[0,10*2*pi/S.freq(min(index(3)))],'figure',gcf); % Plot frequency spectrum of z_1 of specific index
solplot_output  = S.solplot(DYN,solplot_options);



end