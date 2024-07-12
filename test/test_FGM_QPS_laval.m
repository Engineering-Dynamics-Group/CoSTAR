
    load('workspace_init_laval_QPS_FGM');
    
    mu_limit = [1.75,2.5];
    
    eta = mu_limit(1,2);
    Di = 0.2;
    Delta = 1/3;
    e = 0.25;
    d3 = 0.25;
    Fg = 0.3924;
    
    param = {eta,Di,Delta,e,d3,Fg};
    
    cont = 1;
    active_parameter = 1;
    Fcn = @(t,z,param)laval_qp(t,z,param);
    non_auto_freq = @(mu) mu;
    
    
    options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);    %Properties of the System
    options.opt_sol  = costaropts('stability','off','cont','on','auto_freq',auto_freq,'non_auto_freq',non_auto_freq,'sol_type','quasiperiodic','approx_method','fourier-galerkin','act_param',active_parameter);      %Properties of the solution
    options.opt_cont = costaropts('step_control','angle','direction',-1,'pred','tangent','subspace','pseudo-arc','mu_limit',mu_limit,'step_width',0.01,'max_cont_step',1e4);                                                             %Properties for continuation
    options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','int_poincare');   %Properties for approx_method (e.g. Shoot)
    options.opt_init = costaropts('c0',zeros(4,1),'cmatrix',c_max,'smatrix',s_max,'hmatrix',K3);
    
    tic
    [S,DYN] = costar(options);
    toc
    
    %% Test Postprocessing
    benchmark_postprocess_quasiperiodic(DYN,S);