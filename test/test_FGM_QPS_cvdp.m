
    eps_limit = [0.01,1.5];
    
    eps =   0.01;
    alpha = 0.1;
    beta = 1.1;
    param = {eps,alpha,beta};
    

    cont = 1; 
    active_parameter = 1;
    Fcn = @(t,z,param)coupledvdp(t,z,param);
  
    load('workspace_init_cvdp_QPS_FGM');
    
        %Readjust the higher harmonics matrix
    K3 = [0 1 0  3  0 -1 1;...
          0 0 1  0  3  2 2];

    %% Properties
    options.system   = costaropts('order',1,'rhs',Fcn,'param',param,'dim',4);    %Properties of the System
    options.opt_sol  = costaropts('stability','off','cont','on','auto_freq',auto_freq,'sol_type','quasiperiodic','approx_method','fourier-galerkin','act_param',active_parameter);      %Properties of the solution
    options.opt_cont = costaropts('step_control','angle','direction',1,'pred','tangent','subspace','pseudo-arc','mu_limit',eps_limit,'step_width',0.01,'max_cont_step',1e4);                                                             %Properties for continuation
    options.opt_approx_method = costaropts('n_FFT',2^6,'phasecond','int_poincare','n_hh_max',50);   %Properties for approx_method (e.g. Shoot)
    options.opt_init = costaropts('c0',zeros(4,1),'cmatrix',c_max,'smatrix',s_max,'hmatrix',K3);

    
    tic
    [S,DYN] = costar(options);      
    toc


    %% Test Postprocessing
    benchmark_postprocess_quasiperiodic(DYN,S);
