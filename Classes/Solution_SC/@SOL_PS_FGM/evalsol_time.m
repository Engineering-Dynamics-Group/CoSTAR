% Method of SOL_PS_FGM: This method calculates or evaluates a time dependent solution trajectory(!) for the indices defined in options.index
% This method gets called by the solget method of the superclass Solution  
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s_out:   Time solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @t_out:   Array of the time points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 
    

    
function  [s_out,mu,t_out] = evalsol_time(obj,DYN,options)
    
    
    index   =   options.index;
    N       =   numel(index);                       %Number of solutions asked for
    dim = DYN.dim;                                  %dimension of the state space
    mu = obj.mu(1,index);
    counter = 0;

    %Alloquate memory
    s_out = zeros(options.resolution, DYN.dim, N);
    t_out = zeros(options.resolution, 1, N);
    
    for k = index
        counter = counter + 1;
        n_hh = (size(obj.s{1,k},1)/dim-1)/2+1;               %Compute the number of higher harmonics
    
        %This could also be done in the real domain... however: this is how we
        %do it in the residuum equation
        FC = [obj.s{1,k}(1:dim,1);obj.s{1,k}((dim+1):(n_hh)*dim,1)-1i.*obj.s{1,k}(((n_hh)*dim+1):end)];      %Assemble complex Fourrier vector
        FC = reshape(FC,dim,n_hh);
        freq    = obj.freq(1,k);
    
    
        if isfield(options,'interval')          %If an integration interval was supplied by the user... use this interval - else: used the standard interval
            t = linspace(options.interval(1),options.interval(2),options.resolution);
        else
            t = linspace(0,2*pi/freq,options.resolution);                                                     %Integration interval
        end
    
        %Resolution is different from n_FFT and might change... recompute the complex_harmonic_functions:
    
        chf        =  exp(1i.*freq.*obj.hmatrix{1,k}'.*t);        %complex harmonic functions: e.g. real(obj.p_chf(2,:)) gives the cosine with the first higher harmonic defined in obj.hmatrix (first one is 0)
        %imag(obj.p_chf(2,:)) correspondingly gives the same for the sine function
    
        s_out(:,:,counter) = permute(real(pagemtimes(FC,chf)),[2,1,3]);
        t_out(:,1,counter) = t.';
    end
    
end