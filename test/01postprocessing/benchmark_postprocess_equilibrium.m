function benchmark_postprocess_equilibrium(DYN,S)


%% Contplot
contplot_options = costaropts('zaxis','max2','resolution',200,'Color','r');                     % Plot bifurcation diagram using max(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis','mean2','resolution',200,'figure',gcf,'Color','b');       % Plot bifurcation diagram using mean(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis','min2','resolution',200,'figure',gcf,'Color','g');        % Plot bifurcation diagram using min(norm())
contplot_output  = S.contplot(DYN,contplot_options);

contplot_options = costaropts('zaxis',@(z) z(1),'resolution',200,'Color',[0,0,0]);              % Plot bifurcation diagram with user-set function
contplot_output  = S.contplot(DYN,contplot_options);


end