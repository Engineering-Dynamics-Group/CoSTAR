% This function is a method of Solution subclass SOL_PS_FDM and is called by the solget method of the superclass Solution
% It (re-)calculates a FFT of the quasi-periodic orbit data of a solution parametrised in time.
% The resolution of the time analysis is options.resolution. 
% The FFT gives back the absolute value of the complex fourier transform and applies a rectangular leackage window.
% ATTENTION: This method must be adapted when error control for finite differences is implemented
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s:       Amplitude array: This must(!) be a [options.resolution/2 x state_space_dimension x n_evals] dimensional array !!!
% @a:       Phase angle array: This must(!) be a [options.resolution/2 x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @f:       Array of the frequency points: This must(!) be a [options.resolution/2 x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 



function [s,a,mu,omega] = evalsol_frequency(obj,DYN,options)


    [s_time,mu,time] = obj.evalsol_time(DYN,options);                   % Get the time dependent solutions
    

    % Parameters
    dim = DYN.system.dim;                                               % Dimension of state space
    n_evals = length(options.index);                                    % Number of evaluations
    

    % Preallocating f, s and a
    omega = zeros(floor(options.resolution/2), 1, n_evals);             % floor(options.resolution/2) due to output of costarFFT
    s = zeros(floor(options.resolution/2), dim, n_evals);               % floor(options.resolution/2) due to output of costarFFT
    a = zeros(floor(options.resolution/2), dim, n_evals);               % floor(options.resolution/2) due to output of costarFFT


    % Loop for the evaluations
    for i = 1:n_evals

        [f,s_amp,s_angle]   =  costarFFT(time(:,1,i),s_time(:,:,i));        % Compute the Fourier-Transform
        
        omega(:,1,i) = 2*pi*f;
        s(:,:,i) = s_amp;
        a(:,:,i) = s_angle;

    end


end