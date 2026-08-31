function M = expandResidualMatr(obj,input_matrix,dim,n_char,n_shoot)
%EXPANDRESIDUALMATR
% Embeds the characteristic boundary residual into the complete
% multiple-shooting residual matrix.
%
% input_matrix:
%   (dim*n_char) x (dim*n_char+1)
%
% M:
%   (dim*n_shoot*n_char) x (dim*n_shoot*n_char+1)

n_boundary = dim*n_char;
n_state = dim*n_shoot*n_char;

first_column = input_matrix(:,1);
square_matrix = input_matrix(:,2:end);

% M = zeros(n_state,n_state+1,'like',input_matrix);

% Offset of every characteristic group in the complete MS vector.
offsets = (0:n_char-1)*(n_shoot*dim);

% Columns belonging to the first shooting node of each characteristic.
first_node_columns = reshape( ...
    (1:dim).' + offsets,[],1) + 1;

% Remaining shooting-node columns.
other_node_columns = reshape( ...
    (dim+1:n_shoot*dim).' + offsets,[],1) + 1;

% Rows occupied by the boundary residuals.
boundary_rows = reshape( ...
    ((n_shoot-1)*dim + (1:dim)).' + offsets,[],1);

% Insert derivatives with respect to the first shooting nodes.
M(boundary_rows,first_node_columns) = square_matrix;

% All structurally unrelated perturbations retain the base residual.
n_other_columns = n_state-n_boundary;

if n_other_columns > 0
    M(boundary_rows,other_node_columns) = ...
        repmat(first_column,[1,n_other_columns]);
end

% First column contains the unperturbed residual.
M(boundary_rows,1) = first_column;

end