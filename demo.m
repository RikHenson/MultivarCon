close all
clear;
clc;
rng('default')

% Parameters setting
nSubj = 10;
nt = 200; %number of time points
nvoxels = 100; %number of voxels (must be multiple of 2 for below)
nruns = 2; %number of runs. MVPD requires at least 2 runs.
opt.method='pca_exvar';
opt.percentage=99;

%% First example: presence of anticorrelated voxel activities within a ROI

fnam = 'Presence of anticorrelated voxel activities within a ROI';
T = eye(nvoxels); % mapping across ROIs
C = kron([1 -0.99; -0.99 1],eye(nvoxels/2)); % correlation within ROI
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nt,nvoxels),C);
        Yb{g}{r} = Ya{g}{r}*T + randn(nt,nvoxels);
    end
end
plotmv(fnam,T,Ya,Yb,opt,1);

%% Second example: anticorrelation within a ROI induced by the functional mapping
fnam = 'Anticorrelation within a ROI induced by the functional mapping';

T = kron([1 -1; -1 1],eye(nvoxels/2));
C = kron([1 0.1; 0.1 1],eye(nvoxels/2)); % correlation within ROI
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nt,nvoxels),C);
        Yb{g}{r} = Ya{g}{r}*T + 0.8*randn(nt,nvoxels);
    end
end
plotmv(fnam,T,Ya,Yb,opt,1);


%% Third example: presence of structured noise

fnam = 'Presence of structured noise'; 
T = eye(nvoxels);
C = kron([1 0.9; 0.9 1],eye(nvoxels/2));
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nruns
        Ya{g}{r} = mvnrnd(zeros(nt,nvoxels),C);
        Yb{g}{r} = Ya{g}{r}*T + randn(nt,nvoxels);
        Yb{g}{r} = Yb{g}{r}   + 5*repmat(randn(nt,1),1,nvoxels);
    end
end
plotmv(fnam,T,Ya,Yb,opt,0);


%% Fourth example: stimulus-dependent linear mapping

fnam = 'Stimulus-dependent linear mapping';
C = eye(nvoxels);
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:3
       Ya{g}{r} = mvnrnd(zeros(nt,nvoxels),C);
       T = randn(nvoxels);
       Yb{g}{r} = Ya{g}{r}*T + randn(nt,nvoxels);
    end
end
plotmv(fnam,T,Ya,Yb,opt,0);

%% Fifth example: pooling across multiple instances of a condition makes RC more sensitive at the RC level for independent noise

clear Ya Yb
% identical trial-averaged activity patterns
pat_common = randn(1,nvoxels);
ntrials = 20;nruns = 5;
Ya{1} = repmat(pat_common,[ntrials,1]);Yb{1} = repmat(pat_common,[ntrials,1]);
n1 = 5*randn(size(Ya{1}));n1 = n1-mean(n1);
n2 = 5*randn(size(Yb{1}));n2 = n2-mean(n2);
Ya{1} = Ya{1} + n1;Yb{1} = Yb{1} + n2;

% trial-averaged RC (correlation of condition x condition RDMs)
rc_cond = corr(mean(Ya{1})',mean(Yb{1})');
% trial-level RC (correlation of trial x trial RDMs)
rc_trial = corr(pdist(Ya{1})',pdist(Yb{1})');
fprintf('RC at the condition level: %.2f\nRC at the trial level: %.2f\n',rc_cond,rc_trial)


return