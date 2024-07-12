% No function yet!!!!!!
%
%   This is a method of the Stability subclass ST_QPS_SHM
%   It calculates the Lyapunov exponents and indicates if the solution is stable or unstable
%
%@obj:  ST_QPS_SHM object
%@y:    current solution point
%@J:    Jacobian of the current solution point
%@DYN:  Dynamical System class object
%@AM:   Approximation Method subclass object
%
%@multipliers: Lyapunov exponents
%@n_unstable:  number of unstable Lyapunov exponents
%@max_mult:    Lyapunov exponent closest to stability boundary (centered)

function [multipliers,n_unstable,max_mult] = QPS_calc_stability(obj,y,J,DYN,AM)

global SubJacobianSysGlob InputArrayGlob
persistent IdxEV IdxREALout IdxIMAGout

% Variables
LenZ = prod(InputVar.NumOfElem);
LenZX = prod(InputVar.NumOfElem+1);
NOZ = InputVar.NumberOfVariables.NumberOfAmplitudes;
IdxZ = 1:LenZ:NOZ*LenZ;

% Time derivatives of finite differences
DerivZ = SystemOfPDE(InputArrayGlob,InputVar.Para);
NOZpP = size(DerivZ,1);

% Index of jacobian entries w/o phase
IdxJMatrix = reshape(1:NOZpP^2,NOZpP,NOZpP);
IdxJ = reshape(IdxJMatrix(1:NOZ,1:NOZ),1,[]);

% Preallocation
JacobianArrayTemp = zeros(LenZ,NOZpP.^2);
ValueJaciPhase = zeros(NOZpP,InputVar.ProdNumOfElem);
JacobianFunction{NOZpP,NOZpP} = [];
FundMatrixFUN{NOZpP,NOZpP} = [];

%% Calculate Linearized System
% Jacobian at every point along the N x 1 vector.
% Perturbe given phases
InputArrayGlobPer = sqrt(eps).*(1+abs(InputArrayGlob(NOZ+1:end,:)));

% Perturbe all column                                          % Jacobian
% of system for linearized system
for j=1:size(InputArrayGlobPer,1)
    DerivZPer = SystemOfPDE([InputArrayGlob(1:NOZ,:);InputArrayGlob(NOZ+1:end,:)+InputArrayGlobPer(j,:)],InputVar.Para);
    ValueJaciPhase(:,:,j) = 1./InputArrayGlobPer(j,:).*(DerivZPer-DerivZ);
end

% Colums follow the matlab reshape logic
for i=1:LenZ
    ValueJaci = reshape(SubJacobianSysGlob(IdxZ+i-1,IdxZ+i-1),1,[]);
    % All variables w/o non-autonomous frequencies
    for j=1:NOZ^2
        k = IdxJ(j);
        JacobianArrayTemp(i,k) = ValueJaci(j);
    end
    if NOZ ~= NOZpP
        % Complete non-autonomous entries
        JacobianArrayTemp(i,IdxJMatrix(1,NOZ+1):end) = reshape(ValueJaciPhase(:,i,:),1,[]);
    end
end

% Consider periodicity
JacobianMatrixDummy = reshape(JacobianArrayTemp,[InputVar.NumOfElem,NOZpP.^2]);
JacobianMatrix = cat(1,JacobianMatrixDummy,JacobianMatrixDummy(1,:,:));
if size(InputVar.NumOfElem,2)>1
    JacobianMatrix = cat(2,JacobianMatrix,JacobianMatrix(:,1,:));
end
JacobianArray = reshape(JacobianMatrix,[LenZX,NOZpP.^2]);

% Interpolate Data
idx = reshape(1:NOZpP^2,[NOZpP NOZpP]);

% Define function for jacobian
% Preallocation
T1C = 0:(2*pi/InputVar.NumOfElem(1)):2*pi;
T2C = 0:(2*pi/InputVar.NumOfElem(2)):2*pi;
[T2,T1] = meshgrid(T2C,T1C);
ThC = cat(3,T1,T2);

% Define start values for Theta: All coordinates on the first boundary (Number of characteristics)
NOC = InputVar.NumberOFCharacteristics;
% Define start values for DeltaZ: Identity matrix
IC = repmat(reshape(eye(NOZpP),[],1),NOC,1);

% Preallocation
ThetaInArra = zeros(2,NOC);

% Check if frequency of quasiperiodic motion is non-autonomous or autonomous
% (CAUTION: Higher frequency is chosen)
if isempty(IntFreeFreq)
    FreqTheta =  [DerivZ(end-1,1);DerivZ(end,1)];
else
    % Case: One autonomous and one non-autonomous
    if numel(IntFreeFreq) == 1
        FreqTheta =  [DerivZ(end,1);IntFreeFreq(1)];
        % Case: Two autonomous
    else
        FreqTheta =  [IntFreeFreq(1);IntFreeFreq(2)];
    end
end
% The lower frequency is the discretization direction, NOT the higher one !
[FreqValMax,FreqIdxMax] = max(FreqTheta);
[~,FreqIdxMin] = min(FreqTheta);
TimeIntFreq = FreqValMax;
ThetaInArra(FreqIdxMin,:) = linspace(0,2*pi*(NOC-1)/NOC,NOC);

% Jacobian function
for h=1:NOZpP
    for i=1:NOZpP
        JacobianFunction{i,h} = griddedInterpolant(T1,T2,JacobianMatrix(:,:,IdxJMatrix(i,h)),'spline');
    end
end

% Jacobian as matrix function
Jacobian = @(x)cellfun(@(c) c(x),JacobianFunction);

% Calculate Fourier series of Jacobian enties
FourierSeriesFunction = CalcFSQuasiperiodic(Jacobian,JacobianArray,NOZpP,ThC);

% Time Simulation
TimeMax = abs(2*pi/TimeIntFreq);
TimeSpan = [0 TimeMax];
SimOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'MaxStep',1e-1);

% Cell for block reshape
BlockIDx{1,size(ThetaInArra,2)} = [];
for i=1:size(ThetaInArra,2)
    BlockIDx{1,i} = (i-1)*NOZpP+1:i*NOZpP;
end

% Time simulation
disp('Time integration...')
tic
[~,FundMatrixDummy] = ode45(@(t,z) TimeSimStability(t,z,FourierSeriesFunction,NOZpP,FreqTheta,ThetaInArra,BlockIDx),TimeSpan,[IC;zeros(size(InputVar.NumOfElem,2),1)],SimOpts);
toc

FundMatrix = reshape(FundMatrixDummy(end,1:end-InputVar.DimNumOfElem),NOZpP,NOZpP,[]);

% Save data for branch switchting
if InputVar.IndiBranchSwitchingTool
    save(['./../Case_',num2str(InputVar.Count),'/0DATA_NSB_',num2str(ContPara,'%.8f'),'.mat'])

    % Empty (indicator) persitent variables
    IdxEV = [];
end

%% -----------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Post-processor
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Parameter
DataSamp = 10000;
Indi = 1;
p = 64;

% Consider periodicity
FundMatrixARRAY  = cat(3,FundMatrix,FundMatrix(:,:,1));

% Function of the fundamental solution along Theta1
for j=1:NOZpP
    for i=1:NOZpP
        FundMatrixFUN{i,j} = griddedInterpolant([ThetaInArra(FreqIdxMin,:),2*pi],permute(FundMatrixARRAY(i,j,:),[1,3,2]),'spline');
        %                 FundMatrixFUN{i,j} = griddedInterpolant([ThetaInArra(FreqIdxMin,:),2*pi],FundMatrixARRAY(i,j,:),'spline');
    end
end
% Modulate function as array output
FundFUN = @(x)cellfun(@(c) c(x),FundMatrixFUN,'UniformOutput',false);

% Calculate Ljapunov Exponents
% Number of Multiplication
NOM = 2*10^4;
x = 1:NOZpP;

% Take triangular matrix instead of identity matrix for columns in the jacobian with only one entry not equal to zero
FundNorm = eye(NOZpP);

disp('Mapping...')
tic

% All coordinates of mapping
TCR = permute(mod(abs(((1:NOM)-1)*2*pi*(FreqTheta(FreqIdxMin))/(FreqTheta(FreqIdxMax))),2*pi),[1 3 2]);
% All fundamental matrices of mapping
FundARRAY = cell2mat(FundFUN(TCR));
FundARRAYVal = zeros(size(FundARRAY));

for i=1:NOM
    % Multiplication: current fundamental soultion and normalized solution
    FundARRAYVal(:,:,i) = FundARRAY(:,:,i)*FundNorm;

    % Gram-schmidt orthonormalisation
    FundNorm = GSOrthonormalization(FundARRAYVal(:,:,i));
end

TimesFunc = @(x) FundARRAYVal(:,:,x)'*FundARRAYVal(:,:,x);
Vsqr = reshape(cell2mat(arrayfun(TimesFunc, 1:NOM,'UniformOutput',false)),NOZpP,NOZpP,[]);      % Ersetzen

% Calculate volume
detFunc = @(x,y) sqrt(det(Vsqr(1:x,1:x,y)));
Vol = reshape(arrayfun(detFunc,repmat(x,1,NOM),reshape(repmat(1:NOM,NOZpP,1),1,[])),NOZpP,[])';

toc

% Catch errror
if max(max(imag(Vol))) > sqrt(eps)
    error('!!! USERINFO: Time interval is too large for stability algorithm .... !!!')
end
% Lyapunov exponents of n-th order
LyExN = real(sum(log(Vol)))./(NOM*TimeSpan(2));

% Calculate LEs of first order
Output.Values = sort([LyExN(1),diff(LyExN)],'descend');
Output.Info = [];

end




