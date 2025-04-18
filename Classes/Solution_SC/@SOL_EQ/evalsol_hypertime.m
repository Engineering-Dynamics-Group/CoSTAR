%Method of Equilibrium: This method returns the solution vector s as well as a 
%meaningless hypertime array (equal to hypertime array of periodic solutions) for all requested indices
%
%@obj:      obj of clas Equilibrium < Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        torus solution array: 
%           !!! This must(!) be a [solution_dimension x conti_steps] dimensional array !!!
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x conti_steps] dimensional array !!!
%@hypertime:vector of the evaluated hypertime coordinates:
%           !!! This must(!) be a [options.resolution x 1 x conti_steps] dimensional array !!! dimensional array !!!


function [s_out,mu,hypertime] = evalsol_hypertime(obj,DYN,options)

    s_out = obj.s(:,options.index);
    mu  = obj.mu(1,options.index);
    hypertime = [];
    
end