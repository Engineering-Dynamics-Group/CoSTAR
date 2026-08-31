function M = perturbMatrSub(obj,input_matrix,dim,n_char,n_shoot)
%PERTURBMATRSUB
% Arranges the perturbed values of all shooting nodes in a block-diagonal
% structure.
%
% input_matrix:
%   dim x n_shoot x n_char x (dim+1)
%
% M:
%   (dim*n_shoot*n_char) x (dim*n_shoot*n_char+1)

n_blocks = n_shoot*n_char;
n_values = dim*n_blocks;

values = reshape(input_matrix,[n_values,dim+1]);

% Initially, every column contains the unperturbed values.
M = repmat(values(:,1),[1,n_values+1]);

% Because MATLAB linearizes [dim,n_shoot,n_char], the shooting-node
% index changes first and the characteristic index second.
for i_block = 1:n_blocks

    rows = (i_block-1)*dim + (1:dim);
    columns = rows + 1;

    M(rows,columns) = values(rows,2:dim+1);

end

end