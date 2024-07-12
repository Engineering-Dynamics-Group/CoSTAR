% Function defines subspace constraint and derivative for continuation
% algorithm
%
%@obj:      Continuation class object
function obj = choose_subspace(obj)

if strcmpi(obj.subspace,'pseudo-arc')                                       %Pseudo-arc-length continuation
    obj.sub_con = @(y,obj) obj.dy0.'*(y-obj.yp);
    obj.d_sub_con = @(y,obj) obj.dy0.';
elseif strcmpi(obj.subspace,'natural')                                      %Natural step subspace constraint
    obj.sub_con = @(y,obj) y(end)-obj.yp(end);
    obj.d_sub_con = @(y,obj) [zeros(1,size(y,1)-1),1];
elseif strcmpi(obj.subspace,'arclength')                                    %Radial (arclength) subspace constraint
    obj.sub_con = @(y,obj) (y-obj.y0).'*(y-obj.y0)-obj.step_width.^2;
    obj.d_sub_con = @(y,obj) 2.*(y-obj.y0).';
elseif strcmpi(obj.subspace,'taxi')                                         %1-Norm subspace constraint
    obj.sub_con = @(y,obj) norm(y-obj.y0,1)-obj.step_width;
    obj.d_sub_con = @(y,obj) (sign(y-obj.y0)).';
end

end