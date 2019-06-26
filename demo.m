close all
clear;
clc;
rng('default')

% Parameter settings
nSubj = 20;         % number of subjects (replications)
nTime = 200;        % number of time poinTimes
nVoxs = [50 60];    % number of voxels in ROI1 and ROI2
mVoxs = min(nVoxs);
nruns = 2;          % number of runs. MVPD requires at least 2 runs.
Noise = 1;          % std of Noise in ROI2

% MVPD parameters
opt.method='pca_exvar';
opt.percentage=99;

%% First example: positively correlated voxel activities within ROI1 (where both UniConn and MultiConn work)

fnam = 'Positively correlated voxel activities within a ROI';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = [1 2 1 2];
%vis = 2;
plotmv(fnam,T,C,Ya,Yb,opt,vis);

%% Second example: anticorrelated voxel activities within ROI1 (where MulitCon work better)

fnam = 'Negatively correlated voxel activities within a ROI';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.9;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + (1-cc)*eye(nVoxs(1));
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
plotmv(fnam,T,C,Ya,Yb,opt,vis);

%% Third example: anticorrelation in ROI2 induced by the functional mapping (where MulitCon work better)

fnam = 'Negative correlations induced by the functional mapping';
%T = rand(nVoxs)-0.5;
T = zeros(nVoxs); for j=1:mVoxs/2; T(j,j)=1; T(j+mVoxs/2,j+mVoxs/2)=-1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
plotmv(fnam,T,C,Ya,Yb,opt,vis);

%% Fourth example: run-dependent linear mapping (where MVPD fails)

fnam = 'Run-dependenTime linear connectivity';
cc = 0.5;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
%C = eye(nVoxs);
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
       Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       T = rand(nVoxs)-0.5;
       Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
plotmv(fnam,T,C,Ya,Yb,opt,vis);

%% Fifth example: the functional mapping is nonlinear (where Dcor works best)

fnam = 'Nonlinear connectivity';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
       Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       Yb{g}{r} = abs(Ya{g}{r}*T) + Noise*randn(nTime,nVoxs(2));
    end
end
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
vis = 2;
plotmv(fnam,T,C,Ya,Yb,opt,vis);

%% Sixth example: the presence of structured noise in ROI2 (where PCA works best)

fnam = 'Structured Noise in ROI2'; 
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
        Yb{g}{r} = Yb{g}{r}   + 5*Noise*repmat(randn(nTime,1),1,nVoxs(2));
    end
end
vis = 2;
plotmv(fnam,T,C,Ya,Yb,opt,vis);


% %% Fifth example: pooling across multiple instances of a condition makes RC more sensitive at the RC level for independenTime Noise
% 
% clear Ya Yb
% % idenTimeical trial-averaged activity patterns
% pat_common = randn(1,nVoxs);
% nTimerials = 20;nruns = 5;
% Ya{1} = repmat(pat_common,[nTimerials,1]);Yb{1} = repmat(pat_common,[nTimerials,1]);
% n1 = 5*randn(size(Ya{1}));n1 = n1-mean(n1);
% n2 = 5*randn(size(Yb{1}));n2 = n2-mean(n2);
% Ya{1} = Ya{1} + n1;Yb{1} = Yb{1} + n2;
% 
% % trial-averaged RC (correlation of condition x condition RDMs)
% rc_cond = corr(mean(Ya{1})',mean(Yb{1})');
% % trial-level RC (correlation of trial x trial RDMs)
% rc_trial = corr(pdist(Ya{1})',pdist(Yb{1})');
% fprinTimef('RC at the condition level: %.2f\nRC at the trial level: %.2f\n',rc_cond,rc_trial)


return