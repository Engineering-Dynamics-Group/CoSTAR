%This is a general purpose method for sorting the multipliers (eigenvalues,Floquet,Lyapunov) based on different criteria 
%
%@obj:          Stability subclass object
%@DYN:          DynamicalSystem object
%@multipliers   Eigenvalues, Floquet mulitpliers or Lyapunov Exponents
%@vectors       Vectors corresponding to multipliers
%
% AS FOR NOW: The function is inteded to be expanded by more advanced sorting criteria (MAC,...)
function varargout = sort_multipliers(obj,DYN,varargin)

    switch numel(varargin)
        
        case 1
            multipliers = varargin{1};
        
        case 2
            multipliers = varargin{1};
            vectors = varargin{2};
    end


    switch DYN.sol_type

        case 'equilibrium'

            [~,idx] = sort(real(multipliers));      %Sort the eigenvalue according to their real value

        case 'periodic'

            [~,idx] = sort(multipliers);       %Sort the Floquet multipliers according to their absolute value

        case 'quasi-periodic'

            [~,idx] = sort(multipliers);            %Sort the Lyapnuov according to their  value

    end


    %Sort the multipliers, the vectors (if present) and define the output
    switch nargout

        case 1
            varargout{1} = multipliers(idx);

        case 2
            if numel(varargin) ~= 2             %Make sure that multipliers and vectors were passed to this method
                warning('Error using Stability method "sort_multipliers": multipliers and vectors have been defined as output, but there is an input missing.')
            end
            varargout{1} = multipliers(idx);
            varargout{2} = vectors(:,idx);

    end


end