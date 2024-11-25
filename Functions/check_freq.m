% This function checks if provided frequency(s) or frequency(s) of a computed solution are >= freq_limit
%
% @freq_limit:  Bottom limit for the frequency(s)
% @freq: Frequency(s) (scalar or vector)
% @DYN:  DynamicalSystem object
% @y:    Solution vector
% @argin_1: - Call from Gatekeeper: @freq_limit
%           - Call from other method: @DYN
% @argin_2: - Call from Gatekeeper: @freq
%           - Call from other method: @y
% @warn_msg:    Warning message


function warn_msg = check_freq(argin_1,argin_2)

    formatSpec = '%.0e';                        % Set format of number output
    warn_msg = '';                              % Initialize warning message as empty char array


    % This case is primarily meant to be used by the Gatekeeper
    if isa(argin_1,'double')

        freq_limit = argin_1;                   % Frequency limit is the first input argument
        freq       = argin_2;                   % Frequency (vector) is the second input argument

        if any(freq < freq_limit)               % Check if any frequency is smaller than the defined limit
            warn_msg = append('CoSTAR stopped because at least one frequency is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');   % Rather "simple" error message, but it is not used anyway
            % If it is required, the error messages could be defined more precisely (see elseif below).
            % This could be achieved by defining the input variable as struct with fields Omega, omega, Omega_1, omega_1, Omega_2 and omega_2.
            % Depending on which field(s) exist, we could distinguish between the different cases and therefore set a precise error message.
        end



    % This case is used when check_freq is called in initial_solution and in m_continuation
    elseif isa(argin_1,'DynamicalSystem')

        DYN    = argin_1;                       % DynamicalSystem object
        y      = argin_2;                       % Solution vector
        n_freq = DYN.n_freq;                    % Number of frequencies of the solution
        n_auto = DYN.n_auto;                    % Number of autonomous frequencies of the solution
        freq_limit = DYN.freq_limit;            % Get frequency limit from DYN


        % Periodic solutions
        if n_freq == 1

            if n_auto == 0                      % Non-autonomous system
                mu = y(end);                    % Get continuation parameter
                Omega = DYN.non_auto_freq(mu);  % Get non-autonomous frequency
                if Omega < freq_limit
                    warn_msg = append('CoSTAR stopped because the non-autonomous frequency at the predictor point or of the latest solution equals Omega = ', num2str(Omega), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                end

            elseif n_auto == 1                  % Autonomous system
                omega = y(end-1);               % Get autonomous frequency
                if omega < freq_limit
                    warn_msg = append('CoSTAR stopped because the autonomous frequency at the predictor point or of the latest solution equals omega = ', num2str(omega), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                end
            end


        % Quasi-Periodic solutions
        elseif n_freq == 2

            if n_auto == 0                      % Non-autonomous system
                mu = y(end);                    % Get continuation parameter
                Omega = DYN.non_auto_freq(mu);  % Get non-autonomous frequencies
                if (Omega(1) < freq_limit) && (Omega(2) < freq_limit)
                    warn_msg = append('CoSTAR stopped because the non-autonomous frequencies at the predictor point or of the latest solution equal Omega_1 = ', num2str(Omega(1)), ' and Omega_2 = ', num2str(Omega(2)), ', which are below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                elseif Omega(1) < freq_limit
                    warn_msg = append('CoSTAR stopped because the first non-autonomous frequency at the predictor point or of the latest solution equals Omega_1 = ', num2str(Omega(1)), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                elseif Omega(2) < freq_limit
                    warn_msg = append('CoSTAR stopped because the second non-autonomous frequency at the predictor point or of the latest solution equals Omega_2 = ', num2str(Omega(2)), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                end

            elseif n_auto == 1                  % Half-Autonomous (mixed) system
                mu = y(end);                    % Get continuation parameter
                Omega = DYN.non_auto_freq(mu);  % Get non-autonomous frequency
                omega = y(end-1);               % Get autonomous frequency
                if (Omega < freq_limit) && (omega < freq_limit)
                    warn_msg = append('CoSTAR stopped because both frequencies at the predictor point or of the latest solution equal Omega = ', num2str(Omega), ' and omega = ', num2str(omega), ', which are below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                elseif Omega < freq_limit
                    warn_msg = append('CoSTAR stopped because the non-autonomous frequency at the predictor point or of the latest solution equals Omega = ', num2str(Omega), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                elseif omega < freq_limit
                    warn_msg = append('CoSTAR stopped because the autonomous frequency at the predictor point or of the latest solution equals omega = ', num2str(omega), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                end

            elseif n_auto == 2                  % Non-autonomous system
                omega = y(end-2:end-1);         % Get autonomous frequencies
                if (omega(1) < freq_limit) && (omega(2) < freq_limit)
                    warn_msg = append('CoSTAR stopped because the autonomous frequencies at the predictor point or of the latest solution equal omega_1 = ', num2str(omega(1)), ' and omega_2 = ', num2str(omega(2)), ', which are below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                elseif omega(1) < freq_limit
                    warn_msg = append('CoSTAR stopped because the first autonomous frequency at the predictor point or of the latest solution equals omega_1 = ', num2str(omega(1)), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                elseif omega(2) < freq_limit
                    warn_msg = append('CoSTAR stopped because the second autonomous frequency at the predictor point or of the latest solution equals omega_2 = ', num2str(omega(2)), ', which is below the frequency limit of ', num2str(freq_limit,formatSpec), '.');
                end
            end

        end


    end


end