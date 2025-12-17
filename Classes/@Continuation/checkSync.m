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
Z = reshape(s,[DYN.dim,AM.n_int_1,AM.n_int_2]);   % Reshape the solution vector (dim x n_int_1 x n_int_2)

DeltaMu = y1(end,1)-S.mu(1,end);

indFL = indicatorLocking(OMEGA,DeltaMu);

% Check stopping Criterion
if(indFL<1e-3)
    obj.p_contDo = 0;
    obj.p_stopping_flag = 'CoSTAR stopped because solution synchronized by Frequency Locking!';
end


end


% Indicator function for frequency locking
function inFL = indicatorLocking(OMEGA,DeltaMu)

omB(1,:) = OMEGA(1,:)-OMEGA(2,:);
dom_dmu = 1./DeltaMu.*(omB(1,2)-omB(1,1));
inFL = (1./dom_dmu).^2;

end