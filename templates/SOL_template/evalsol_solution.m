%Method of SOL_template: This method simply extracts the solution quantities
%of the SOL_template class object. 
%
%@obj:      obj of clas SOL_template < Solution
%@DYN:      DynamicalSystem object
%@options:  options structure for postprocessing solutions
%
%@s:        time solution container: 
%           !!! This must(!) have [solution method dimension x solution_dimension x
%           conti_steps]  !!! This must (in general) not be an array. This can also be a
%           cell, ...
%@mu:       vector of the evaluated continuation parameters:
%           !!! This must(!) be a [1 x
%           conti_steps] dimensional array !!!



function [s,mu,options] = evalsol_solution(obj,DYN,options)

            %Example:
%             s(:,:,:) = permute(obj.s(:,options.index),[3,1,2]);    %options.index is unique due to S.solget_up_index
%             mu(1,:)  = obj.mu(1,options.index);     %options.index is unique due to S.solget_up_index

  

        
end