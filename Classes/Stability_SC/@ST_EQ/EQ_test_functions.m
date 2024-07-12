%This function evaluates  test_functions for determination of the type of a local co-dimension 1 bifurcation 
%and is used as well for the determination of its position. The test function is assigned to the
%ST.test_function handle in the s_stability_selection method
%
%@obj:          object of Stability subclass EQ_ST_PS
%@multipliers:  Eigenvalues of current solution point
%@DYN:          DynamicalSystem class object
%
%@tfcn:     values of the test function for the current solution point.


function       tfcn = EQ_test_functions(obj,multipliers,DYN)             
   
    %The product of diag with the upper triangular matrix produces the products of the Floquet-Multipliers
    num = numel(multipliers)-1;
  
    A = triu(repmat(fliplr(multipliers(1:num).'),num,1)) + triu(repmat(flipud(multipliers(2:num+1)),1,num));
   
    tmp = [prod(multipliers);                                    %Test for fold bifurcation
    prod(A(triu(true(size(A)))))];                               %Testing for Hopf bifurcation: triu(...) gives a logical matrix for extracting only the upper triangular elements which we want in the product    
    if max(abs(imag(tmp)))>1e-14; warning('Imaginary part detected within PS_test_functions. Contact a CoSTAR programmer.'); end
    tfcn = real(tmp);


end