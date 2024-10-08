% Method of SOL_PS_SHM: This method (re-)calculates a periodic orbit of a shooting solution
% starting from the iterated point on the periodic orbit
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
  
    [s_hypertime,mu,time] = obj.evalsol_time(DYN,options);

    % Normalizes time to [0, 2*pi]
    hypertime = time.*permute(obj.freq(1,options.index),[1,3,2]);
    
end
