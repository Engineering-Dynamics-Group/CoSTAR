%Method of SOL_template: This method calculates the hypertime manifolds of
%a solution, which is 0-/1-/2- or n-dimensional
%
%
%@obj:      obj of clas SOL_template < Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        hypertime solution array: 
%           !!! This must(!) be a [options.resolution x options.resolution  x solution_dimension x
%           conti_steps] dimensional array !!!
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x
%           conti_steps] dimensional array !!!
%@hypertime:        hypertime array: 
%           !!! This must(!) be a [options.resolution x options.resolution x
%           conti_steps] dimensional array !!!

function [s_hypertime,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

  
        

    
    
end
