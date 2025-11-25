%Function is a method of Approximation Method subclass AM_PS_FGM and computes
%an initial value vector base on the opt_init option structure.
% If there are less initial conditions than higher harmonics in obj.hmatrix,
% additional intial values are guessed. There are no real limitations on the initial ordering of obj.hmatrix.
% the only requirement is, that the initial conditions for the Fourier coefficients correspond to the positions
% in the hmatrix, which also might contain more frequencies (which must not necessarily be "higher").
%
% Example: C0 = zeros(2,1), Cmatrix = Smatrix = [1,2;3,4], hmatrix = [0,1,3, 7,2,5,11,9].
%
%@obj:  Object of AM_PS_FGM
%@DYN:  DynamicalSystem class object

function   obj = getIV(obj,DYN)

    %Either one of the two is present (ensured by the gatekeeper)
    if isempty(obj.fc0)
        obj.fc0 = [obj.c0(:);obj.cmatrix(:);obj.smatrix(:)];
    end
    
    [obj.iv,obj.hmatrix] = obj.sort_guess_FC(DYN,obj.fc0);                  %This functionality gets also used by the IF_increase_discretization method
    
    if ~isempty(DYN.auto_freq)                                              %system is (partly) autonomous
        obj.iv = [obj.iv;DYN.auto_freq(:)];                                 %Append the initial_condition vector, if it is an autonomous system.
    end

end