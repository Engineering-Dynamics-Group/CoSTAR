%Method of Equilibrium: This methods "calculates" the FFT of an Equilibrium

%@obj:      obj of class PS_Shoot_Sol < Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        frequency solution array: 
%           !!! This must(!) be a [options.resolution x solution_dimension x conti_steps] dimensional array !!!
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x conti_steps] dimensional array !!!
%@f:        array of the frequency points:
%           !!! This must(!) be a [options.resolution x 1 x conti_steps]  dimensional array !!!

function [s,mu,f] = evalsol_frequency(obj,DYN,options)

       [s_time,mu,time] = obj.evalsol_time(DYN,options);

       %Compute the Fourier-Transform
       for k = 1:numel(options.index)
           t_step = (max(time(:,1,k))-min(time(:,1,k)))./(length(time(:,1,k))-1);
           Fs = 1/t_step;
           L = size(time,1); 

           f_tmp = Fs/L*(0:(L/2-1));

           f(:,1,k) = 2.*pi.*f_tmp.';
           s(:,:,k) = zeros(numel(f_tmp),size(s_time,2));
           s(1,:,k) = s_time(1,:,k);
       end

end