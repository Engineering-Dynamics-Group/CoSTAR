%This is a method of superclass Stability.
%It adds a new curve point of a bifurcation point iteration and its stability information into the curve container
%at the appropriate position and deletes the appropriate enty of the container according to its size p_container_size
%This acts as an advanced set function of the container with FIFO principle.
%The entry gets deleted in such a way that the "new" obj.curve_container allows an identification of the zero-crossing of
% max_multiplier
%
%The arclength of the newly added point is computed in this function. 
%
%@obj:                  Stability subclass object
%@max_multiplier:       maximal multiplier responsibel for stability property of solution
%@y:                    curve point
%@n_unstable:           number of unstable multipliers
%
% TODO: This code needs to be cleaned up and simplified - however, it works this way
function update_curve_container_bfp(obj,y,multipliers,n_unstable)

    n_unstable_arr = cell2mat(obj.curve_container(3,:));          %Find index of stability change 
    idx = find(diff(n_unstable_arr));                             %New point needs to be inserted after index idx

    tmp = sign(diff(cell2mat(obj.curve_container(2,:))));       %directions of the curve for computing the arclength into the right direction
    if ~(all(tmp == 1)||all(tmp == -1)); error('Curve direction changes.'); end %Short safety break for developers: Ensures that all signs are the same.
    direction = tmp(1);             %Guess the CON.direction parameter of the curve from data   
    
    arcl_1 = obj.curve_container{2,idx}+direction.*norm(obj.curve_container{1,idx}-y);  %Compute arclength of point to be added

    tmp2 = cell(5,size(obj.curve_container,2)+1);   %Can't use obj.p_container_size, since the container might not be fully filled  

    tmp2(1,:)     = {obj.curve_container{1,1:idx},y                                 ,obj.curve_container{1,idx+1:end}};
    tmp2(2,:)     = {obj.curve_container{2,1:idx},arcl_1                            ,obj.curve_container{2,idx+1:end}};
    tmp2(3,:)     = {obj.curve_container{3,1:idx},n_unstable                        ,obj.curve_container{3,idx+1:end}};   
    tmp2(4,:)     = {obj.curve_container{4,1:idx},multipliers                       ,obj.curve_container{4,idx+1:end}};
    tmp2(5,:)     = {obj.curve_container{5,1:idx},obj.test_functions(multipliers)   ,obj.curve_container{5,idx+1:end}};

    %Identify the idx of the element in curve_container, which must be deleted:
    %The element before and after the change of stability are excluded and stay
    %in the container (obviously...). We delete the element with the largest index distance
    %to the element before the change of stability.
    idx2     = find(diff(cell2mat(tmp2(3,:))));   
    idx_val  = setdiff(1:size(obj.curve_container,2)+1,[idx2,idx2+1]); % +1 since the newly iterated curve point has been found
    [~,idx3] = max(abs(idx_val-idx2)); %idx3 is the index of the element farest away from the posc. This element is deleted.
    
    tmp2(:,idx_val(idx3)) = [];
    obj.curve_container = tmp2;
   
end