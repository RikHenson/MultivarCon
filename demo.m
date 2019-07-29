close all
clear;
clc;
rng('default')

% Parameter settings
nSubj = 20;         % number of subjects (replications)
nTime = 200;        % number of time points
nVoxs = [50 60];    % number of voxels in ROI1 and ROI2
mVoxs = min(nVoxs);
nRuns = 2;          % number of runs. MVPD requires at least 2 runs.
Noise = 1;          % std of Noise in ROI2

% MVPD parameters
opt.method='svd_exvar';
opt.percentage=90;
opt.number=1;
opt.regularisation=10.^(-1:0.2:3); % regularisation parameter for ridge regression.

%% First example: positively correlated voxel activities within ROI1 (where both UV-conn and MV-conn metrics work)
fnam = 'Positively correlated voxel activities within a ROI';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Y{s}{r} = X{s}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = [1 2 1 2];
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example1.png'),'png')

%% Second example: anticorrelated voxel activities within ROI1 (where MV-conn metrics work better)
fnam = 'Negatively correlated voxel activities within a ROI';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.9;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + (1-cc)*eye(nVoxs(1));
X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Y{s}{r} = X{s}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example2.png'),'png')

%% Third example: anticorrelation in ROI2 induced by the functional mapping (where MVPD and dCor work better)
fnam = 'Negative correlations induced by the functional mapping';
%T = rand(nVoxs)-0.5;
T = zeros(nVoxs); for j=1:mVoxs/2; T(j,j)=1; T(j+mVoxs/2,j+mVoxs/2)=-1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Y{s}{r} = X{s}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example3.png'),'png')

%% Fourth example: run-dependent linear mapping (within-run measures work)
fnam = 'Run-dependent linear connectivity';
cc = 0.5;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
C = cov(rand(nVoxs(1))); %this example doesn't require a particular cov
% structure in either ROI.
X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
       X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       T = randn(nVoxs);% changed this from rand(nVoxs)-.5 
       Y{s}{r} = X{s}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example4.png'),'png')

%% Fifth example: the functional mapping is nonlinear (dCor works best)
fnam = 'Nonlinear connectivity';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
%C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
X = {}; Y = {};
Noise = 0.5;
for s=1:nSubj
    for r=1:nRuns
       X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       Y{s}{r} = abs(X{s}{r}*T) + Noise*randn(nTime,nVoxs(2));
%        Y{s}{r} = (X{s}{r}*T).^2 + Noise*randn(nTime,nVoxs(2));% also a good example
    end
end
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
vis = 2;
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example5.png'),'png')

%% Sixth example: the presence of structured noise in ROI2 (where RCA works best)
fnam = 'Structured Noise in ROI2'; 
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
% T = randn(nVoxs);
cc = 0.5;
%C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
% C = cov(rand(nVoxs(1)));% this can also be general. again we don't need a 
% specfic C for this
X = {}; Y = {};
Noise = 1;
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Y{s}{r} = X{s}{r}*T + Noise*randn(nTime,nVoxs(2));
        Y{s}{r} = Y{s}{r}   + 5*Noise*repmat(randn(nTime,1),1,nVoxs(2));
    end
end
vis = 2;
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example6.png'),'png')

%% Seventh example: averaging timepoints (trials) with same stimulus improves RCA
fnam = 'Negative correlations induced by the functional mapping; averaging across stimuli of same type';
%T = rand(nVoxs)-0.5;
T = zeros(nVoxs); for j=1:mVoxs/2; T(j,j)=1; T(j+mVoxs/2,j+mVoxs/2)=-1; end
Noise = 1;
nStim = 20;
nRep  = nTime/nStim;  % Assumes a factor of nTime
stimuli = repmat([1:nStim],1,nRep);

cc = 0;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI

X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nStim,nVoxs(1)),C);
        X{s}{r} = repmat(X{s}{r},nRep,1);  % Repeat same pattern 
        Y{s}{r} = X{s}{r}*T + Noise*randn(nTime,nVoxs(2));
%        X{g}{r} = X{g}{r} + Noise*randn(nTime,nVoxs(1));  % independent noise
    end
end

vis = 2;
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
% plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
% saveas(gcf,'Graphics/mvcon_example7a.png','png')
rc_orig = mean(MVconn.RCA)-mean(MVconn_null.RCA);

sX = {}; sY = {};
for s=1:nSubj
    for r=1:nRuns
        for stim = 1:nStim
            sX{s}{r}(stim,:) = mean(X{s}{r}(find(stimuli==stim),:));
            sY{s}{r}(stim,:) = mean(Y{s}{r}(find(stimuli==stim),:));
        end
    end
end

vis = 2;
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[sMVconn,sMVconn_null] = computeMVconn(sX,sY,opt);
% plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
% saveas(gcf,'Graphics/mvcon_example7b.png','png')
rc_pooled = mean(sMVconn.RCA)-mean(sMVconn_null.RCA);
fprintf('RC from the original trial x trial RDMs: \t: %.2f\n',rc_orig)
fprintf('RC from the pooled data i.e. stim x stim RDMs: \t: %.2f\n',rc_pooled)

%% Eighth example: closed-loop (where lagged MV-Conn metric works better)
fnam = 'Closed-loop between subpopulations';

% Parameter settings
nTime = 15360;      % number of time points, e.g. 15360 is equivalent to 60 seconds 
                    % acquisition if the resolution is 256 Hz
nVoxs = [12 10];    % number of locations/voxels in ROI1 and ROI2

opt.segleng=256;    % length (in bins) of the segment used for the cross-spectrum computation,
                    % e.g. 256 is equivalent to 1 second segment with a resolution of 256 Hz
opt.freqbins=31:71; % low gamma band, i.e. from 30 Hz to 70 Hz

delaYin= 10;       % delay in bins of the interaction between the regions,
                    % e.g. if the resolution is 256 Hz and the delaYin is 10, 
                    % the actual delay between a and b is delaYin/(256 Hz)= 40ms 
Noise=0.1;

% functional mapping between subpopulations (from ROI1 to ROI2)                 
Tab=randn(floor(nVoxs/2));
% functional mapping between subpopulations (from ROI2 to ROI1)  
Tba=randn(floor(nVoxs(2:-1:1)/2));

% covariance matrix for the subpopulation in ROI1 that leads the one in ROI2
cc = 0.5;
Ca = ones(floor(nVoxs(1)/2))*cc + eye(floor(nVoxs(1)/2))*(1-cc);
% covariance matrix for the subpopulation in ROI2 that leads the one in ROI1
cc = 0.5;
Cb = ones(floor(nVoxs(2)/2))*cc + eye(floor(nVoxs(2)/2))*(1-cc);

% correlated noise due to source-leakage/volume-conduction/field-spread
Mixnoise=randn(nTime,sum(nVoxs));

X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        
        % lagged multivariate interaction from a population in ROI1
        % (called LPa, i.e. leading population in ROI1) to a population in ROI2
        % (called FPb, i.e. following population in ROI2)
        LPa=mvnrnd(zeros(nTime+delaYin,floor(nVoxs(1)/2)),Ca);
        FPb=(LPa(delaYin+1:end,:)*Tab)+Noise*Mixnoise*randn(sum(nVoxs),floor(nVoxs(2)/2));
        LPa=LPa(1:nTime,:)+Noise*Mixnoise*randn(sum(nVoxs),floor(nVoxs(1)/2));
        X{s}{r}=LPa;
        Y{s}{r}=FPb;

        % lagged multivariate interaction from a population in ROI2
        % (called LPb, i.e. leading population in ROI2) to a population in ROI1
        % (called FPa, i.e. following population in ROI1)
        LPb=mvnrnd(zeros(nTime+delaYin,floor(nVoxs(2)/2)),Cb);
        FPa=(LPb(delaYin+1:end,:)*Tba)+Noise*Mixnoise*randn(sum(nVoxs),floor(nVoxs(1)/2));
        LPb=LPb(1:nTime,:)+Noise*Mixnoise*randn(sum(nVoxs),floor(nVoxs(2)/2));
        X{s}{r}=[X{s}{r}, FPa];
        Y{s}{r}=[Y{s}{r}, LPb];

    end
end

vis = [1 2 1 2];
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,Tab,Ca,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example8.png'),'png')

return

