function M = perturbMatr(obj,input_matrix, dim, n_char)
%PERTURBMATR
% Arranges the perturbed characteristic values in a block-diagonal
% structure.
%
% input_matrix:
%   dim x n_char x (dim+1)
%
% M:
%   (dim*n_char) x (dim*n_char+1)
%
% The first column contains the unperturbed values. Each following
% column represents the perturbation of one scalar characteristic state.
n_values = dim*n_char;

values = reshape(input_matrix,[n_values,dim+1]);

% Initially, every column contains the unperturbed values.
M = repmat(values(:,1),[1,n_values+1]);

% Insert the perturbed dim x dim blocks.
for i_char = 1:n_char

    rows = (i_char-1)*dim + (1:dim);
    columns = rows + 1;

    M(rows,columns) = values(rows,2:dim+1);

end

end