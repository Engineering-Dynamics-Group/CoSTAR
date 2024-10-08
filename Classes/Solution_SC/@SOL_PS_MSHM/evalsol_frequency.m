% Method of SOL_PS_SHM: This method (re-)calculates a FFT of the periodic orbit data of a shooting solution
% starting from the iterated point on the periodic oribt. The resolution of
% the time analysis is options.resolution. The FFT gives back the absolute
% value of the complex fourier transform and a applies a rectangular
% leackage window
%
% @obj:     Solution subclass object
% @DYN:     DynamicalSystem object
% @options: options structure for postprocessing solutions
%
% @s:       Amplitude array: This must(!) be a [options.resolution x state_space_dimension x n_evals] dimensional array !!!
% @mu:      Vector of the evaluated continuation parameters: This must(!) be a [1 x n_evals] dimensional array !!!
% @f:       Array of the frequency points: This must(!) be a [options.resolution x 1 x n_evals]  dimensional array !!!
% n_evals:  Number of curve points to be evaluated 

function [s,mu,f] = evalsol_frequency(obj,DYN,options)

    [s_time,mu,time] = obj.evalsol_time(DYN,options);
    
    %Compute the Fourier-Transform
    for k = 1:numel(mu)
        [f_tmp,s_tmp]   =  costarFFT(time(:,1,k),s_time(:,:,k));
        f(:,1,k) = 2.*pi.*f_tmp.';
        s(:,:,k) = s_tmp;
    end

end