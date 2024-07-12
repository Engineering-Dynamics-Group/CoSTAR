%This is a function of superclass Stability.
%If a change of stability happens on a curve, this function 
%approximates this Point Of Stability Change (posc) based 
%on the data stored in the curve_container property. 
%Interpolation of data is done w.r.t. the arclength of the curve.
%
%@obj:  stability subclass object
%@DYN:  DynamicalSystem object
%
%@y:    approximated curve point of stability change

function y_app = approx_posc(obj,DYN)

    y_app = 0;                                              % ???

    stab_fcn = obj.get_stability_fcn(DYN);
    
    
    % Interpolation step: linear interpolation (duh)
    if size(obj.curve_container,2) == 2                     % ???

        y0 = obj.curve_container{1,1};
        y1 = obj.curve_container{1,2};

        y_app = 1/(stab_fcn(2)-stab_fcn(1)) .* (stab_fcn(2).*y0 - stab_fcn(1).*y1);

        % The method below is the more extensive way of calculating y_app - inserting everthing into one equation leads to the formula above
        % arcl = cell2mat(obj.curve_container(2,:));
        % s = diff(fliplr(arcl).*stab_fcn)./diff(stab_fcn);
        % y_app = 1./diff(arcl).*((s-arcl(1)).*y1 - (s-arcl(2)).*y0);
        
    end

end