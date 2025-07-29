%%Function gives back a string, which identifies the detected bifurcation
%
%@obj:          Stability subclass object
%
%@label:        label identifying the bifurcation
%@msg:          message containing more information on the bifurcation found



function [label,msg] = identify_bifurcation(obj,J)

    curve_container = obj.curve_container;
    
    if ~any(isnan(cell2mat(curve_container(5,:))))        %Are there test functions? If not - use the multipliers
        test_fcn     = cell2mat(curve_container(5,:));
        idx     = find(diff(sign(test_fcn),1,2)); 
        
        if ~isempty(idx)    %Is there a sign change in the test functions? If not - use the multipliers
            idx = idx(1);
            idx = obj.AdditionalConstraints(idx,J);
            label = obj.bifurc_label{1,idx};
            msg = obj.msg_label{1,idx};
        else
            label = obj.bifurc_label{1,end};
            msg = obj.msg_label{1,end};

        end
    else
            label = obj.bifurc_label{1,end};
            msg = obj.msg_label{1,end};
    end

end

