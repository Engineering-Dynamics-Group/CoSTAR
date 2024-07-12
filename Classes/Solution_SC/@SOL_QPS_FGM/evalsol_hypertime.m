% Method of SOL_QPS_FGM: This method calculates the hypertime manifolds of solution, which is 2-dimensional
%
% @obj:       Solution subclass object
% @DYN:       DynamicalSystem object
% @options:   options structure for postprocessing solutions
%
% @s_hypertime: Hypertime solution array: This must(!) be a [options.resolution x options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:          Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @hypertime:   Array of the time points: This must(!) be a [options.resolution x options.resolution x 2 x n_evals] dimensional array !!!
% n_evals:      Number of curve points to be evaluated 

function [s_hypertime,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

    index = options.index;
    n_idx = numel(index);                       %Number of solutions asked for
    
    dim = DYN.dim;                          %dimension of the state space
    counter = 0;
 
    %%% THIS LINE MUST BE CHANGED IF HMATRIX IS ADAPTED IN FUTURE RELEASE
    
    %Alloquate memory
    s_hypertime = zeros(options.resolution,options.resolution,dim,n_idx);

    for k = index

        counter = counter + 1;

        s = obj.s{1,k};                     %Matrix of Fourier coefficients      
        hmatrix = obj.hmatrix{1,k};
        n_hh = (size(s,1)/dim-1)/2+1;           %Compute the number of higher harmonics


        %This could also be done in the real domain... however: this is how we
        %do it in the residuum equation
        FC = [s(1:dim,1);s((dim+1):(n_hh)*dim,1)-1i.*s(((n_hh)*dim+1):end,1)];        %Assemble complex Fourrier vector    
        FCtemp = reshape(FC,dim,n_hh);

        %Resolution is different from n_FFT might change... recompute the complex_harmonic_functions:

        tmp = linspace(0, 2*pi, options.resolution);        %evaluation theta_values
            
        hypertime  = [reshape(repmat(tmp' ,[1,options.resolution]),[1,options.resolution^2]);                       
                     reshape(repmat(tmp   ,[options.resolution,1]),[1,options.resolution^2])]; 
            
        %How does this work?: I am evaluating here the arguments in the cosine or sine functions: cos(H1*theta1+H2*theta2). hmatrix'*p_freq_val is the scalar product
        %The definition of p_freq_val (similar to the meshgrid function) ensures, that p_chf is evaluated at every point of the discretized torus. 

        chf        = exp(1i.*hmatrix'*hypertime);        %complex harmonic functions: e.g. real(obj.p_chf(2,:)) gives the cosine with the first higher harmonic defined in obj.hmatrix (first one is 0)
     
                                                                    %imag(obj.p_chf(2,:)) correspondingly gives the same for the sine func
                                                               
   
        s_hypertime(:,:,:,counter) = reshape(real(FCtemp*chf).',[options.resolution,options.resolution,dim]);    

    end

    mu = obj.mu(1,index);
    hypertime = repmat(reshape(hypertime.',[options.resolution,options.resolution,2]),[1,1,1,n_idx]);
    
end
