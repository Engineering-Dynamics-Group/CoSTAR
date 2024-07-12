% Method of SOL_PS_FGM: This method calculates the hypertime manifolds of a solution, which is 1-dimensional
%
% @obj:       Solution subclass object
% @DYN:       DynamicalSystem object
% @options:   options structure for postprocessing solutions
%
% @s_hypertime: Hypertime solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:          Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @hypertime:   Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:      Number of curve points to be evaluated 

function [s_hypertime,mu,hypertime] = evalsol_hypertime(obj,DYN,options)


    index = options.index;
    N = numel(index);                       %Number of solutions asked for
    counter = 0;

    hypertime = linspace(0, 2*pi, options.resolution);      %evaluation theta_values
    mu = obj.mu(1,index);

    dim = DYN.dim;                          %dimension of the state space
   
    %Alloquate memory
    s_hypertime = zeros(options.resolution,DYN.dim,N);



    for k = index


        counter = counter + 1;
        n_hh = (size(obj.s{1,k},1)/dim-1)/2+1;                                               %Compute the number of higher harmonics
        
        %This could also be done in the real domain... however: this is how we
        %do it in the residuum equation
        
        FC = [obj.s{1,k}(1:dim,1);obj.s{1,k}((dim+1):(n_hh)*dim,1)-1i.*obj.s{1,k}(((n_hh)*dim+1):end)];      %Assemble complex Fourrier vector    
        FC = reshape(FC,dim,n_hh);
    
        %Resolution is different from n_FFT. hmatrix might change for every steo
        chf        = exp(1i.*obj.hmatrix{1,k}'*hypertime);              %complex harmonic functions: e.g. real(obj.p_chf(2,:)) gives the cosine with the first higher harmonic defined in obj.hmatrix (first one is 0)  
                                                                                    %imag(obj.p_chf(2,:)) correspondingly gives the same for the sine function 
        s_hypertime(:,:,counter) = permute(real(pagemtimes(FC,chf)),[2,1,3]);    
        

    end
    hypertime = repmat(hypertime.',[1,1,N]);
     
end
