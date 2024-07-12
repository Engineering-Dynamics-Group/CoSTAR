%This is a method of the Stability superclass and gets used by all subclasses.
%It checks whether one of the input arguments contains inf or NaN values. If 
%true it returns all values set to NaN to indicate the user, that something went 
%wrong.
%
%obj:               Stability subclass object
%multipliers:       Eigenvalues, Floquet or Lyapunov values
%vectors:           Vectors corresponding to multipliers
%n_unstable:        indicates the number of unstable exponents


function  [multipliers,vectors,n_unstable] = check_stability_values(obj,multipliers,vectors,n_unstable,stabi_flag)


        if any(isinf(multipliers))||any(isnan(multipliers))||(stabi_flag==0)
            multipliers = NaN(size(multipliers));
            vectors = NaN(numel(multipliers));
            n_unstable = NaN;
        end

end