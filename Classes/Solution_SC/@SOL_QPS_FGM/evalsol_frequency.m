% Method of SOL_QPS_FGM: This method does evaluate the Fourier spectrum directly
% The method returns the absolute value and the phase angle of the complex fourier transform 
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @P:       Amplitude array: This must(!) be a [options.resolution/2 x state_space_dimension x n_evals] dimensional array !!!
% @a:       Phase angle array: This must(!) be a [options.resolution/2 x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @f:       Array of the frequency points: This must(!) be a [options.resolution/2 x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function [P,a,mu,f] = evalsol_frequency(obj,DYN,options)

    index   = options.index;
    N       = numel(index);         % Number of solutions asked for
    dim     = DYN.dim;              % Dimension of the state space
    res     = options.resolution;   % Resolution 
    counter = 0;

    %Initialise:
    P = zeros(res/2,dim,N);
    a = zeros(res/2,dim,N);
    f = zeros(res/2,1,N);
    
    % Task: define a frequency vector with the demanded resolution, which includes all known frequencies freq.*hmatrix! 
    % In other words: the frequencies apart from freq.*hamtrix are filler frequencies with an amplitude of 0.
    
    % TODO: Vectorise this code for better performance.

    for k = index

        counter = counter + 1;
        s = obj.s{1,k};
        n_hh = (size(s,1)/dim-1)/2+1;               % Compute the number of higher harmonics
  

        % This could also be done in the real domain... however: this is how we do it in the residuum equation
        FC          = [s(1:dim,:);s((dim+1):(n_hh)*dim,:)-1i.*s(((n_hh)*dim+1):end,:)];      % Assemble complex Fourrier vector    
        FC          = reshape(FC,dim,n_hh);
        freq        = obj.freq(:,k);
        hmatrix     = obj.hmatrix{1,k};             % Harmonics

        [frq_sort,idx_sort] = sort(abs(freq.'*hmatrix));            % Some frequency combination might be negative. frq_sort thus must be sorted in ascending order
        A = abs(FC(:,idx_sort));                                    % Order the fourier-coefficients accordingly to the frequencies and compute the absolute amplitudes
        alpha = atan2(-imag(FC(:,idx_sort)),real(FC(:,idx_sort)));  % Order the fourier-coefficients accordingly to the frequencies and compute the phase angles

        if res > n_hh

            f_val = sort([frq_sort(2:end),linspace(0,1.1*max(frq_sort),res/2-numel(frq_sort(2:end)))]);   % Add as many zeros as are needed for the the vector to have the length resolution  
            f_val_u = unique(f_val);                % Remove the doubled frequencies (in most cases, there should be not more than one since it is very unlikely that linspace hits a value of frq_sort(2:end-1))
            % Most of the values of f_val_u are equidistantly distributed within [0,max(frq_sort)]. However, the frequencies ...
            % frq_sort are placed somewhere in between, which can lead to ugly plots. That is why the frequency values between ...
            % the frq_values (overall numel(frq_values)-1 sections) are distributed equidistantly within each section again.
            [~,idx] = intersect(f_val_u,[frq_sort,1.1*max(frq_sort)]);          % Find the index positions where the frequencies frq_sort are located
            for i = 1:(numel(idx)-1)
                f_val_u(idx(i):idx(i+1)) = linspace(f_val_u(idx(i)),f_val_u(idx(i+1)),numel(idx(i):idx(i+1)));
            end
            % Now that the frequency values are distributed equidistantly as far as possible, we need to add some frequencies in case that numel(f_val_u) ~= res/2
            spacing = f_val_u(end) - f_val_u(end-1);        % Compute a mean spacing for adding frequencies
            % spacing = mean(diff(f_val));                  % OLD CODE: spacing
            f(:,1,counter) = [f_val_u, max(f_val_u) + (1:(numel(f_val)-numel(f_val_u)))*spacing].';     % Add the number of removed frequencies at the end with the appropriate spacing to restore the correct length

            % Now add the amplitude and angle values at the appropiate positions. 
            % We can use idx(1:end-1) for that (f(idx(end)) = 1.1*max(frq_sort))
            % [~,idx] = intersect(f(:,1,counter),frq_sort); % OLD CODE: Find the indices of the contained frequencies    
            P(idx(1:end-1),:,counter) = A.';                % Write the values
            a(idx(1:end-1),:,counter) = alpha.';            % Write the values
            
        else

            f(:,1,counter) = freq.*hmatrix(1:res);
            P(:,:,counter) = A(:,1:res).';
            warning(append('Request of frequency for Fourier-Galerkin method: Your requested options.resolution = ',num2str(res),...
                           ' is smaller than the number of frequencies = ', num2str(numel(hmatrix)),'. Consider raising the resolution.'));
        
        end

    end

    mu(1,:)  = obj.mu(1,index);  

end