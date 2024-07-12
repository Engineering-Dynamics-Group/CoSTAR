%This is a method of superclass Stability
%After iterating a bifurcation point, it deletes all elements after 
%the stability change in the container. That way, no "old" point of 
%stability change is contained in the container
%
%@ obj: object of Stability Subclass

function clean_curve_container(obj)
    
        idx = find(diff(cell2mat(obj.curve_container(3,:))));                       %Idx of element in front of stability change
        
        if ~isempty(idx)
            obj.curve_container(:,1:idx) = [];                                      %Delete all elements in front of point of stability change
        else
            obj.curve_container = cell(0);                                          %Clear everything: This should actually never be the case                                      
        end

end