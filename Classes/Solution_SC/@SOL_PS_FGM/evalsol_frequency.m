% Method of SOL_PS_FGM: This method synthesizes a spectrum from the Fourier coefficients and DOES NOT USE AN FFT
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

    % Initialise:
    P = zeros(res/2,dim,N);
    a = zeros(res/2,dim,N);
    f = zeros(res/2,1,N);

    % Task: define a frequency vector with the demanded resolution, which includes all known frequencies freq.*hmatrix! 
    % In other words: the frequencies apart from freq.*hamtrix are filler frequencies with an amplitude of 0.
   
    for k = index

        counter = counter + 1;
        s = obj.s{1,k};
        n_hh = (size(s,1)/dim-1)/2+1;               % Compute the number of higher harmonics
   
        % This could also be done in the real domain... however: this is how we do it in the residuum equation
        FC          = [s(1:dim,1);s((dim+1):(n_hh)*dim,1)-1i.*s(((n_hh)*dim+1):end,1)];      % Assemble complex Fourrier vector    
        FC          = reshape(FC,dim,n_hh);
        freq        = obj.freq(1,k);
        hmatrix     = obj.hmatrix{1,k};             % Harmonics

        A = abs(FC);                                % Absolute amplitudes
        alpha = atan2(-imag(FC),real(FC));          % Phase angles
                
        if res > n_hh

            % The idea is to depict the frequency domain in the range of [0 freq*(max(hmatrix)+1)].
            % To ensure that the frequency domain contains exactly the values freq.*hmatrix of the amplitude peaks AND has an ...
            % equidistant distribution of the frequency values, we need to find the highest number <= (res/2-1) that is divisible ...
            % by (max(hmatrix)+1) without remainder. The variable n_overshoot tells the difference between (res/2-1) and that number.
            % (The "-1" is to exclude the value at freq*(max(hmatrix)+1) from the search, since we want to distribute ...
            % (res/2-1-n_overshoot) points in (max(hmatrix)+1) "intervals".)
            % After that, we can create the frequency values in the range of [0 freq*(max(hmatrix)+1)] using (res/2-n_overshoot) points.
            % (Be aware that we must not include the "-1" now because we need the value at freq*(max(hmatrix)+1) in the frequency domain.)
            % Finally, we add n_overshoot points for frequency values > (max(hmatrix)+1) using the spacing of f_val to ensure that ...
            % the  frequency domain has exactly res/2 elements.
            n_overshoot = mod(res/2-1,max(hmatrix)+1);
            f_val = linspace(0,freq*(max(hmatrix)+1),res/2-n_overshoot);            % Create (most of) the frequency values
            spacing = f_val(end)-f_val(end-1);                                      % Spacing of f_val
            f(:,1,counter) = [f_val, max(f_val) + (1:n_overshoot)*spacing];         % Add n_overshoot frequency values to ensure that size(f,1) = res/2
            % OLD CODE: 
            % It has the problem that the frequency domain is computed by distributing res/2-numel(hmatrix) evaluation points equidistantly ...
            % in the interval [0,freq.*max(hmatrix)] and placing the numel(hmatrix) harmonic frequencies in between. This can lead to an ...
            % inconsistent distribution of the frequency values near the freq*hmatrix, which results in ugly plots.
            %{
            f_val = sort([freq.*hmatrix,linspace(0,freq.*max(hmatrix),res/2-numel(hmatrix))]);          % Add as many zeros as are needed for the the vector to have the length resolution
            spacing = mean(diff(f_val));            % Compute a mean spacing
            f_val_u = unique(f_val);                % Remove the doubled frequencies
            f(:,1,counter) = [f_val_u, max(f_val_u) + (1:(numel(f_val)-numel(f_val_u)))*spacing].';     % Add the number of removed frequencies at the end with the appropriate spacing to restore the correct length
            %}

            % Now add the amplitude and angle values at the appropiate positions: Find the index positions first
            idx_spacing = (res/2-1-n_overshoot) / (max(hmatrix)+1);                 % Index spacing between two harmonics
            idx = 1 + idx_spacing.*hmatrix;                                         % Thanks to the code above, we know exactly where the harmonics are located within f
            % Alternative code: use the find() function to validate idx
            %{
            idx_val = zeros(1,numel(hmatrix));
            for i = 1:numel(hmatrix)
                idx_val(i) = find(abs(f(:,1,counter)-freq*hmatrix(i)) < sqrt(eps));     % We need to use this due to floating point errors of linspace
            end
            % [~,idx] = intersect(f(:,1,counter),freq.*hmatrix);        % OLD CODE: Find the indices of the contained frequencies (can't be used anymore due to floating point errors of linspace)
            %}
            P(idx,:,counter) = A.';                                     % Write the values
            a(idx,:,counter) = alpha.';                                 % Write the values

        else

            f(:,1,counter) = freq.*hmatrix(1:res);
            P(:,:,counter) = A(:,1:res).';
            a(:,:,counter) = alpha(:,1:res).';
            warning(append('Request of frequency for Fourier-Galerkin method: Your requested options.resolution = ',num2str(res),...
                            ' is smaller than the number of frequencies = ', num2str(numel(hmatrix)),'. Consider raising the resolution.'));
        
        end

    end

    mu(1,:)  = obj.mu(1,index);                     % options.index is unique due to S.solget_up_index

end