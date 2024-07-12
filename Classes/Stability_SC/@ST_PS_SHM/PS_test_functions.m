%This is a collection function for bifurcation test functions for periodic solutions.
%This is a static method since the ST_PS_Hill method should be able to call this function also.
%The test function is assigned to the ST.test_function handle in the s_stability_selection method
%
%@multipliers:  eigenvalue, Floquet multiplier or Lyapunov exponents
%@tfcn:         values of the various test functions
%@DYN:          DynamicalSystem class
%
% test-functions according to "Marx, Vogt - Dynamische Systeme - Theorie und Numerik - S. 345"
function tfcn = PS_test_functions(multipliers,DYN)                  


    if DYN.n_auto == 1

        stab_indicator = abs(multipliers);
        [~,idx] = sort(abs(stab_indicator-1),'ascend');              %Sorting of the absolute (values-1) this way always puts the stab_indicators closest to 1 at the first position.
        multipliers(idx(1)) = [];                                    %Deleting the first of the sorted stab_indicators always eliminates the one closest to 1, which is either the one associated to the autonomous frequencies or one close to it.

    end
    
    %The product of diag with the upper triangular matrix produces the products of the Floquet-Multipliers
    num = numel(multipliers)-1;
    A = diag(multipliers(1:num))*triu(repmat(multipliers(2:num+1).',num,1))-triu(ones(num,num));    %Matrix for testing Neimark-Sacker bifurcation


    tmp = [
    prod(multipliers-1);                                         %Testing for Fold-, Pitchfork- or transcritical bifurcation
    prod(multipliers+1);                                         %Testing for Flip (period doubling) bifurcation    
    prod(A(triu(true(size(A)))))];                               %Testing for Neimark-Sacker bifurcation: triu(...) gives a logical matrix for extracting only the upper triangular elements which we want in the product    

    % if max(abs(imag(tmp)))>1e-14; warning('Imaginary part detected within PS_test_functions. Contact a CoSTAR programmer.'); end
    if max(abs(imag(tmp)))>1e-14 
        warning('Imaginary part detected within PS_test_functions. Reducing the step size or increasing the discretization might solve the issue.'); 
    end
    tfcn = real(tmp);
    
end