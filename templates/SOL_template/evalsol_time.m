%Method of SOL_template: This method calculates or evaluates a time dependent solution 
%trajectory(!) for the indices defined in options.index

%@obj:      obj of clas SOL_template < Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        time solution array: 
%           !!! This must(!) be a [options.resolution x solution_dimension x
%           conti_steps] dimensional array !!!
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x
%           conti_steps] dimensional array !!!
%@t:        array of the time points:
%           !!! This must(!) be a [options.resolution x 1 x
%           conti_steps]  dimensional array !!!

%This method gets called by the solget method of the superclass Solution

function  [s,mu,t] = evalsol_time(obj,DYN,options)
    
    %Define the indices at which solution shall be evaluated:
    
    
    %Define s: Example for Equilibrium
    %s(:,:,:) = permute(obj.s(:,options.index),[3,1,2]);
    
    %Define t
    

    %Get the mu values

    %Example:
    %mu = obj.mu(options.index);    %options.index is unique due to S.solget_up_index



end