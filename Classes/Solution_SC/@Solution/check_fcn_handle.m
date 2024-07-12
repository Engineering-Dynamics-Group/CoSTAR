%This method is a a small submethod for the SOLUTION class gatekeepers checks the function
%handles
%
%@obj:          SOLUTION subclass object
%@GC:           GATEKEEPER class object
%@fcn_handle:   function_handle to be checked
%@name_handle:  name 'string' of the function_handle for error message
%@d_in:         dimension of input of the fcn_handle: can either be the
%strings 'solution_argument' 'matrix','column_vector','row_vector','scalar' (of dimension is not known)
%or an array specifying the dimension. 'solution_argument' automatically selects the right dimension for an argument 
% of the corresponding solution type (eq,p,qp): 1D for eq, 2D for periodic, 3D, for quasi-periodic.
%@d_out:        excpected dimension of the output of the function handle: can either be the
%strings 'matrix','column_vector','row_vector','scalar' (if dimension is not known)
%or an array specifying the dimension 
%
%

function check_fcn_handle(obj,GC,fcn_handle,name_handle,d_in,d_out)


%Determine the solution type (equilibrium, periodic, 2D-quasi-periodic):
if ~isprop(obj,'freq')
    n_freq = 0;
else
    n_freq = size(obj.freq,1);
end


     if isa(d_in,'char')
         if strcmpi(d_in,'solution_argument')
             if n_freq == 0 %equilibrium
                 d_in = [10,1];  
             elseif n_freq == 1 %periodic
                 d_in = [10,10];
             else %quasiperiodic
                d_in = [10,10,5];
             end

         elseif strcmpi(d_in,'matrix')
             d_in = [10,10];
         elseif strcmpi(d_in,'column_vector')
             d_in = [10,1]; %10 is a random number
         elseif strcmpi(d_in,'row_vector')
            d_in = [1,10]; %10 is a random number
         elseif strcmpi(d_in,'scalar')
            d_in = [1,1]; 
         else
            error('Something went wrong');
         end
     end

     if isa(d_out,'char')

         if strcmpi(d_out,'solution_argument')
             if n_freq == 0 %equilibrium
                 d_out = [1,1];  
             elseif n_freq == 1 %periodic
                 d_out = [10,1];
             else %quasiperiodic
                d_out = [10,10];
             end

         elseif strcmpi(d_out,'matrix')
             d_out = [10,10];
         elseif strcmpi(d_out,'column_vector')
             d_out = [10,1]; %10 is a random number
         elseif strcmpi(d_out,'row_vector')
            d_out = [1,10]; %10 is a random number
         elseif strcmpi(d_out,'scalar')
            d_out = [1,1]; 
         else
            error('Something went wrong');
         end
     end

     if ~(nargin(fcn_handle)==1)    %has the function_handle the correct number of input argument
         GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' has ',num2str(nargin(fcn_handle)),' arguments, but is only allowed to have 1.');
     else %function handle has the correct number of input arguments
            
        try
            tmp = fcn_handle(ones(d_in));
            d_fcn = size(tmp);

            if ~(numel(d_fcn)==numel(d_out))    
                    GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns for an input which is a ', num2str(numel(dtmp)),'-dimensional array. However, the output should be a ', num2str(numel(d_out)),'-dimensional array.'); 
            elseif ~all(d_fcn==d_out)
                    GC.error_msg{1,end+1} = append('Your provided function handle ',name_handle,' returns for an example-input of size [', num2str(d_in(1)),',',num2str(d_in(2)),'] an output of size [', num2str(d_fcn(1)),',',num2str(d_fcn(2)),']. However, the output should be of the size [', num2str(d_out(1)),',',num2str(d_out(2)),'].'); 
            end
            
        catch
                    GC.error_msg{1,end+1} = append('Something with your provided function handle ',name_handle,' is completely wrong. It should have the following stats: E.g.: Input size: [', num2str(d_in(1)),',',num2str(d_in(2)),']. --> Output size: [', num2str(d_out(1)),',',num2str(d_out(2)),']');
        end


     end



end