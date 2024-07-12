%Method of Equilibrium: This method gives back a non-sense time dependent evaluation of the equilibrium
%
%@obj:      obj of clas Equilibrium < Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        time solution array: 
%           !!! This must(!) be a [resolution x solution_dimension x conti_steps] dimensional array !!!
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x conti_steps] dimensional array !!!
%@t:        array of the time points:
%           !!! This must(!) be a [resolution x 1 x conti_steps]  dimensional array !!!

function  [s_out,mu,t] = evalsol_time(obj,DYN,options)

        s = obj.s(:,options.index);
        mu(1,:)  = obj.mu(1,options.index);

        if isfield(options,'interval')          %If an integration interval was supplied by the user... use this interval - else: used the standard interval
            tspan = linspace(options.interval(1),options.interval(2),options.resolution);
        else
            tspan = linspace(0,2.*pi,options.resolution);                                                     %Integration interval
        end
        s_out = repmat(permute(s,[3,1,2]),options.resolution,1,1);
        t = repmat(tspan.',1,1,numel(options.index));    %Give back a meaningless time structure in the correct dimension. 

end