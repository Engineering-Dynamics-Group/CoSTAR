%Function is a method of ApproxMethod subclass AM_QPS_FGM and computes
%an initial value vector base on the opt_init option structure
% If there are less initial conditions than higher harmonics in obj.hmatrix,
% additional intial values are guessed
%
%@obj:  Object of AM_QPS_FGM
%@DYN:  DynamicalSystem class object
    
function   obj = getIV(obj,DYN)
   
    %Either one of the two (fc0 or C0/cmatrix/smatrix) is present (ensured by the gatekeeper)
    if isempty(obj.fc0)
        obj.fc0 = [obj.c0(:);obj.cmatrix(:);obj.smatrix(:)];
    end
    
    [obj.iv,obj.hmatrix] = obj.sort_guess_FC(DYN,obj.fc0);                  %This functionality gets also used by the IF_increase_discretization method
    
    if ~isempty(DYN.auto_freq)                                              %system is (partly) autonomous
        obj.iv = [obj.iv;DYN.auto_freq(:)];                                 %Append the initial_condition vector, if it is an autonomous system.
    end

end