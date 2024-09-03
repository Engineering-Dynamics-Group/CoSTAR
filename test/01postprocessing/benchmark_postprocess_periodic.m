function benchmark_postprocess_periodic(DYN,S)

n = size(S.mu,2);                                                                                   % Get number of curve points
index = [round(n/4,0),round(n/2,0),round(3*n/4,0)];                                                 % Get three indecies to plot different solutions


%% Solget
solget_options = costaropts('eval',@(z) z(:,1),'space','time','resolution',1e3,'index',1);
solget_output  = S.solget(DYN,solget_options);

solget_options = costaropts('eval',@(z) max(z),'space','hypertime','index',index);
solget_output  = S.solget(DYN,solget_options);



%% Contplot
contplot_options = costaropts('zaxis','max2','resolution',200,'Color','r');                             % Plot bifurcation diagram using max(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis','mean2','resolution',200,'figure',gcf,'Color','b');               % Plot bifurcation diagram using mean(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis','min2','resolution',200,'figure',gcf,'Color','g');                % Plot bifurcation diagram using min(norm())
contplot_output  = S.contplot(DYN,contplot_options);


contplot_options = costaropts('zaxis',@(z)max(z(:,1)),'resolution',200,'Color',[0,0,0]);                % Plot bifurcation diagram with user-set function
contplot_output  = S.contplot(DYN,contplot_options);


%{
contplot_options = costaropts('zaxis','max2','resolution',200);
contplot_output  = S.contplot(DYN,contplot_options);
%}



%% Solplot
solplot_options = costaropts('zaxis',@(z) z(:,1),'space','hypertime','resolution',200,'index',index(3),'Color','r');                % Plot solution z_1 in hypertime
solplot_output  = S.solplot(DYN,solplot_options);

%{
solplot_options = costaropts('zaxis','all','space','frequency','resolution',100,'index',index(1));
solplot_output  = S.solplot(DYN,solplot_options);
%}

solplot_options = costaropts('zaxis',@(z) z(:,1),'space','time','resolution',100,'index',index);                                    % Plot time solution z_1
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,1),'space','time','index',index(2),'interval',[0,10],'linestyle','--','figure',gcf);  % Plot interval of time solution z_1 of specific index
solplot_output  = S.solplot(DYN,solplot_options);


solplot_options = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','index',index);                                                                       % Plot trajectory of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','mu',S.mu(1,index(2)),'figure',gcf);                                                  % Plot trajectory of specific mu-value
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z)z(:,2),'xaxis',@(z)z(:,1),'space','trajectory','index',index(2),'Color','r');                                                          % Plot trajectory of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,2),'xaxis',@(z) z(:,1),'space','trajectory','resolution',104,'index',index(2),'interval',[0,2],'linestyle','--','figure',gcf);    % Plot interval of trajectory of specific index
solplot_output  = S.solplot(DYN,solplot_options);


solplot_options = costaropts('zaxis','all','space','frequency','index',index(2),'resolution',2^10,'interval',[0,10*2*pi/S.freq(index(2))]);                     % Plot frequency spectrum of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,1),'space','frequency','index',index(1),'resolution',2^10,'interval',[0,10*2*pi/S.freq(index(1))]);               % Plot frequency spectrum of z_1 of specific index
solplot_output  = S.solplot(DYN,solplot_options);

solplot_options = costaropts('zaxis',@(z) z(:,1),'space','frequency','index',index(3),'resolution',2^10,'interval',[0,10*2*pi/S.freq(index(3))],'figure',gcf);  % Plot frequency spectrum of z_1 of specific index
solplot_output  = S.solplot(DYN,solplot_options);



end