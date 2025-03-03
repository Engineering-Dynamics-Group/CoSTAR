% This is a method of the Stability subclass ST_QPS_SHM
% It calculates the derivative dZ~_0/dZ_0 which is used to calculate the
% mapping matrix of quasi-periodic shooting solutions from the Jacobian
%
% @obj:      ST_QPS_SHM object
% @y:        current solution point
% @Omega:    Frequencies of current solution point
% @AM:       Approximation Method subclass object
% @index:    Index of quasi-periodic shooting method to determine integration time
%
% @J:        Matrix dZ~_0/dZ_0

function J = jacobi_int(obj,y,Omega,AM,index)

phi = AM.phi;
dim = AM.n;
n_char = AM.n_char;
Ik = [0,2.*pi./Omega(1,index(1,1))];                                                        % Set time span for integration
Theta = mod(phi(index(1,2),:) + Ik(2)*Omega(1,index(1,2)),2*pi);                            % Calculate end values of characteristics and map back to 0,2pi square
[B0(1,:),~] = sort(Theta);                                                                  % Sort remapped values in ascending order

Z0 = reshape(y(1:dim*n_char),dim,n_char);
dx = sqrt(eps)*(1 + max(abs(y(1:end-1,1)),[],1));                                           % Calculate differential to calculate Jacobian
INIT = repmat(Z0,[1,1,dim+1]);                                                              % Set initial values
for kk = 1:dim
    INIT(kk,:,kk+1) = Z0(kk,:) + dx;                                                        % Set initial values of perturbated characteristics
end

F = repmat(reshape(INIT(:,:,1),[dim*n_char,1]),[1,dim*n_char+1]);                           % Set F to unperturbated initial values
kjj = 0;
for lp = 1:dim:dim*n_char
    kjj = kjj+1;   
    F(lp:lp+dim-1,lp+1:lp+dim) = permute(INIT(:,kjj,2:dim+1),[1,3,2]);                      % Set F to values of perturbated initial values on initial interval [0,2pi(1-1/nchar)]
end

G = reshape(F,[dim,n_char,dim*n_char+1]);

Z0_int = zeros(size(G));
for k=1:dim*n_char+1
    INT = csape([phi(index(1,2),:),2*pi],[G(:,:,k),G(:,1,k)],'periodic');                   % Interpolate the perturbated values of Z
    Z0_int(:,:,k) = fnval(INT,B0(1,:));                                                     % Evaluate the values on the nodes
    clear INT;
end

G_int = reshape(Z0_int,[dim*n_char,dim*n_char+1]);                                          % Reshape interpolated values

J = (G_int(:,2:end)-G_int(:,1))./dx;                                                        % Calculate numeric differential

end