% This method saves or reset properties of the Continuation and ApproxMethod object that are modified by the error control
% If the error control fails to recompute a solution with increased or decreased discretisation, the modified properties ...
% need to be reset in order to comply with the last computed solution. If they are not reset, MATLAB errors concerning array dimensions can occur
%
% @obj:            Continuation class object
% @AM:             ApproxMethod class object
% @save_or_reset:  char determining whether to save or to reset the properties

function obj = ec_save_reset_props(obj,AM,save_or_reset)

    fieldnames_AM_struct = fieldnames(AM.ec_prop_save);
    
    switch save_or_reset

        case 'save'

            % First save the properties of the Continuation object which are modified by the error control (we could also code this ...
            % ... dynamically like below, but since there are only three properties to be saved/reset, it is clearer to do it this way)
            obj.p_ec_prop_save.yp = obj.yp;
            obj.p_ec_prop_save.y0 = obj.y0;
            obj.p_ec_prop_save.dy0 = obj.dy0;

            % Now dynamically save the properties of AM that are modified by the error control
            for i = 1:numel(fieldnames_AM_struct)
                AM.ec_prop_save.(fieldnames_AM_struct{i}) = AM.(fieldnames_AM_struct{i});
            end


        case 'reset'

            % First reset the properties of the Continuation object which are modified by the error control (we could also code this ...
            % ... dynamically like below, but since there are only three properties to be saved/reset, it is clearer to do it this way)
            obj.yp = obj.p_ec_prop_save.yp;
            obj.y0 = obj.p_ec_prop_save.y0;
            obj.dy0 = obj.p_ec_prop_save.dy0;

            % Now dynamically reset the properties of AM that were modified by the error control
            for i = 1:numel(fieldnames_AM_struct)
                AM.(fieldnames_AM_struct{i}) = AM.ec_prop_save.(fieldnames_AM_struct{i});
            end
            
    end

end