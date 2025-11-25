%This is a method of superclass Stability.
%It adds a new curve point and stability information into the curve container
%and deletes the oldest entry of the container according to its size p_container_size
%This acts as an advanced set function of the container with FIFO principle
%
%@obj:          Stability subclass object
%@DYN:          Dynamical system object
%@AM:           Approximation subclass object
%@arcl_1:       arclength of curve point y 
%@y:            curve point
%@mulitpliers:  eigenvalues, Floquet multipliers or Lyapunov multipliers
%@n_unstable:   number of unstable multipliers

function update_curve_container(obj,DYN,AM,arcl_1,y,y0,multipliers,n_unstable)

    %This update works always, even if curve_container is empty   
    obj.curve_container{1,end+1} = y;
    obj.curve_container{2,end}   = arcl_1;
    obj.curve_container{3,end}   = n_unstable;
    obj.curve_container{4,end}   = multipliers;
    obj.curve_container{5,end}   = obj.test_functions(multipliers);

    if size(obj.curve_container,2) > obj.p_container_maxsize %FIFO update of the curve_container
        obj.curve_container(:,1) = [];      %Currently, it should be sufficient to delete only the first solution in the curve container because it is appended by exactly one solution only here in this method (above)
    end

    %If the error_control mechanism is on, we cannot guarantee that all solution curve points in obj.curve_container have the same size (which is necessary for the following interpolation step)
    if strcmpi(AM.error_control,'on') == 1    

        cc_dim = size(obj.curve_container,2);

        if (cc_dim == 2) && (numel(y0) == numel(y))         
            obj.curve_container{1,1} = y0;      %If everything works as intended, y0 and y have the same size. The numel(y0) == numel(y) above is only for safety reasons    
        else
            % Currently, the loops below are unnecessary since curve_container (c_c) only has size(c_c,2) = p_container_maxsize = 2
            % In this case, we could also use the last solution (before y = CON.p_y1) CON.y0 directly since its size was already updated in error_control
            % If p_container_maxsize is changed, the loop becomes necessary again
            for k = 1:size(obj.curve_container,2); cc_dim(1,k) = numel(obj.curve_container{1,k}); end %Get the dimensions of all the solution curve point
            if ~all(cc_dim == cc_dim(1))        %See if all dimensions in the array are the same
                new_dim= cc_dim(end);           %We update the dimension of all curve points in curve_container to the newest one. This is necessary for the method to work with the error_control methods: Specifically: The AM_Q/PS_FGM methods save the last h_matrix used in an error loop, but not all previous h_matrices.
                for idx = 1:(numel(cc_dim)-1)   %Loop through all dimensions smaller than new_dim and update the dimension.
                    obj.curve_container{1,idx} = AM.IF_update_sol_dim(DYN,new_dim,obj.curve_container{1,idx});
                end
            end
        end

    end

end