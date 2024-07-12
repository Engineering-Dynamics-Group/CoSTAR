% This function sets the frequency and intial conditions vectors correctly.
% The function works for equilibrium, periodic and quasi-periodic.
%
%@obj: DynamicalSystem object

function obj = updatefreq(obj)

% obj.ic = obj.opt_init.ic(:);
if  numel(obj.auto_freq) == 2                                   %system is (partly) autonomous
    obj.n_auto = 2;                                             %Set the auto key (determines if system is autonomous/ easier for further computations)
elseif numel(obj.auto_freq) == 1                                %system is (partly) autonomous
    obj.n_auto = 1;  
else
    obj.n_auto = 0;
end

if ~isempty(obj.non_auto_freq)
    if ~isa(obj.non_auto_freq,'function_handle')              %If n_auto_freq is not already a function handle. Make it one.
        obj.non_auto_freq = @(mu) obj.non_auto_freq;
    end
    %Determines the number of base frequencies (this should work for any case since non/auto_freq is allocated with [])
    obj.n_freq = length(obj.auto_freq) + length(obj.non_auto_freq(1)); %We need to call the non_auto_freq handle with a value to determine its output dimension 
else %There is definitely no non-autonomous frequency
    obj.n_freq = length(obj.auto_freq); 
end

end
