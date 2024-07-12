%Method of SOL_template: This method (re-)calculates a FFT of the periodic orbit data of a hypertime solution
%The resolution of the time analysis is options.resolution. The FFT gives back the absolute
%value of the complex fourier transform and a applies a rectangular leackage window
%
%@obj:      obj of clas SOL_template< Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        frequency solution array: 
%           !!! This must(!) be a [options.resolution x solution_dimension x
%           conti_steps] dimensional array !!!
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x
%           conti_steps] dimensional array !!!
%@f:        array of the frequency points:
%           !!! This must(!) be a [options.resolution x 1 x
%           conti_steps]  dimensional array !!!



function [s,mu,f] = evalsol_frequency(obj,DYN,options)

    %Get the torus 
%   [s_hypertime,mu,hypertime] = evalsol_hypertime(obj,DYN,options);
    
    %Compute the Fourier-Transform
    % Do this either by applying the myFFT algorithm or get it directly (Fourier-Galerkin methods)


end