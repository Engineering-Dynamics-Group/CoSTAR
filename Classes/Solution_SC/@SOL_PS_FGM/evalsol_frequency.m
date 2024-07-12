% Method of SOL_PS_FGM: This method synthesizes a spectrum from the Fourier coefficients and DOES NOT USE AN FFT
% The method gives back the absolute value of the complex fourier transform 
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @P:       Frequency solution array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @f:       Array of the frequency points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function [P,mu,f] = evalsol_frequency(obj,DYN,options)

    index   =   options.index;
    N       =   numel(index);                       %Number of solutions asked for
    counter = 0;
    dim = DYN.dim;                          %dimension of the state space

    %Initialise:
    P = zeros(options.resolution,dim,N);
    f = zeros(options.resolution,1,N);

    %Task: define a frequency vector with the demanded resolution, which
    %includes all known frequencies freq.*hmatrix! In other words: the
    %frequencies apart from freq.*hamtrix are filler frequencies with an
    %amplitude of 0.
    
    for k = index

        counter = counter + 1;
        n_hh = (size(obj.s{1,k},1)/dim-1)/2+1;           %Compute the number of higher harmonics
   
    %This could also be done in the real domain... however: this is how we
    %do it in the residuum equation

        FC          = [obj.s{1,k}(1:dim,1);obj.s{1,k}((dim+1):(n_hh)*dim,1)-1i.*obj.s{1,k}(((n_hh)*dim+1):end,1)];        %Assemble complex Fourrier vector    
        FC          = reshape(FC,dim,n_hh);
        freq        = obj.freq(1,k);

        FC = abs(FC);
        hmatrixtmp = obj.hmatrix{1,k};
        hmatrix = obj.hmatrix{1,k};


        if options.resolution > n_hh


            f_val = sort([freq.*hmatrixtmp,linspace(0,freq.*max(hmatrixtmp),options.resolution-numel(hmatrixtmp))]); %Add as many zeros as are needed for the the vector to have the length resolution
            spacing = mean(diff(f_val));    %Compute a mean spacing
            f_val_u = unique(f_val);        %Get out the doubled frequencies
            f(:,1,counter) = [f_val_u, max(f_val_u) + (1:(numel(f_val)-numel(f_val_u)))*spacing].'; %Simply adding the needed number of frequencies at the end with the appropriate spacing in between

            %Now add the FCtmp values at the apropriate positions:
            [~,idx] = intersect(f(:,1,counter),freq.*hmatrix);    %Find the indices of the contained frequencies
            P(idx,:,counter) = FC.';                           %Write the value
        else
            f(:,1,counter) = freq.*hmatrixtmp(1:options.resolution);
            P(:,:,counter) = FC(:,1:options.resolution).';
            warning(append('Request of frequency for Fourier-Galerkin method: Your requested options.resolution = ',num2str(options.resolution),' is smaller than the number of frequencies = ', num2str(numel(hmatrixtmp)),'. Consider raising the resolution.'));
        end
    end

    mu(1,:)  = obj.mu(1,options.index);                    %options.index is unique due to S.solget_up_index

end