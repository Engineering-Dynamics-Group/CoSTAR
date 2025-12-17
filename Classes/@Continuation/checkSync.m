function obj = checkSync(obj,DYN,AM,S)

n_auto = DYN.n_auto;
y1 = obj.p_y1;

if(n_auto==0)
    return;
elseif(n_auto==1)
    OMEGA = [S.freq,[DYN.non_auto_freq(y1(end,1));y1(end-1,1)]];
elseif(n_auto==2)
    OMEGA = [S.freq,y1(end-2:end-1,1)];
end

s = y1(1:(end-(DYN.n_auto+1)),1);
Z = permute(reshape(s,[DYN.dim,AM.n_int_1,AM.n_int_2]),[2,3,1]);            % Reshape the solution vector (dim x n_int_1 x n_int_2)

DeltaMu = y1(end,1)-S.mu(1,end);

indFL = indicatorLocking(OMEGA,DeltaMu);                                    % Evaluate indicator function for locking
indFS_A = indicatorSuppArea(Z);

% Check stopping criterion for Frequency Locking
if(indFL<1e-3)
    obj.p_contDo = 0;
    obj.p_stopping_flag = 'CoSTAR stopped because solution synchronized by Frequency Locking!';
end

% Check stopping criterion for Suppressive Synchronization
[chi_S, idx] = min(indFS_A,[],1);
if(chi_S<1e-1)
    obj.p_contDo = 0;
    obj.p_stopping_flag = append('CoSTAR stopped because solution synchronized by Suppresive Synchronization of Omega_', num2str(idx),'!');
end

% Save values of indicator functions to object of Continuation class
obj.p_IndF_Lock = indFL;
obj.p_IndF_Supp = indFS_A;

end


% Indicator function for frequency locking
function inFL = indicatorLocking(OMEGA,DeltaMu)

omB(1,:) = OMEGA(1,:)-OMEGA(2,:);
dom_dmu = 1./DeltaMu.*(omB(1,2)-omB(1,1));
inFL = (1./dom_dmu).^2;

end


function chi = indicatorSuppArea(Z)

NN = size(Z);
N = NN(1:2);
dim = NN(3);
p = 2;

Zcentered = zeros(N(1,2),N(1,1),dim,p);

c1 = squeeze(sum(Z, 1)) / N(1,1);                                           % Calculate centroids
c2 = squeeze(sum(Z, 2)) / N(1,2);

for l=1:N(1,1)                                                              % Calculate radii with respect to centroids
    Zcentered(:,l,:,1) = permute(Z(l,:,:),[2,3,1,4])-c2(l,:);               % Adjust Poincare sections to be centered by their centroid
end
for l=1:N(1,2)
    Zcentered(l,:,:,2) = permute(Z(:,l,:),[1,3,2,4])-c1(l,:);
end

rCentered(:,:,:) = vecnorm(Zcentered,2,3);                                  % Determine the radius

rAveraged1(:) = mean(rCentered(:,:,1),1);                                   % Average the radii over the firs torus coordinate
rAveraged2(:) = mean(rCentered(:,:,2),2);                                   % Average the radii over the second torus coordinate

rMax1 = max(rAveraged1(:,:),[],2);                                          % Take the maximum of the radii with respect to the other torus coordinate
rMax2 = max(rAveraged2(:,:),[],2);                                          % Take the maximum of the radii with respect to the other torus coordinate

chi = pi.*[rMax2;rMax1].^2;                                                 % Calculate indicator function as areas of max averaged radii of Poincare sections
                                                                            % rMax1 and rMax2 are switched, so that the correct frequency that synchronizes
                                                                            % is identified
end