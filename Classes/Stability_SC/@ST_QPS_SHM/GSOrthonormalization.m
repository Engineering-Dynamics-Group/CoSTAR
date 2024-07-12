%   This is a method of the Stability subclass ST_QPS_SHM
%   It calculates the Gram-Schmidt orthonormalization of a given set of
%   vectors
%
%@obj:      ST_QPS_SHM object
%@MatIn:    Matrix of vectors to be orthonormalized [v1,v2,v3,...] 
%
%@Output:   Matrix with orthonormalized vectors [u1,u2,u3,...]

function Output = GSOrthonormalization(obj,MatIn)
% Gram-Schmidt Orthonormalization
n = size(MatIn,1);
k = size(MatIn,2);
U = zeros(n,k);
            
U(:,1) = MatIn(:,1)/sqrt(MatIn(:,1)'*MatIn(:,1));
for o = 2:k
    U(:,o) = MatIn(:,o);
    for j = 1:o-1
        U(:,o) = U(:,o) - ( U(:,o)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
    end
    U(:,o) = U(:,o)/sqrt(U(:,o)'*U(:,o));
end
Output = U;
end