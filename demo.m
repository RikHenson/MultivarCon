% A Basti, H Nili & R Henson

close all
clear;
clc;

% Parameter settings
nSubj = 20;         % number of subjects (replications)
nTime = 200;        % number of time points
nVoxs = [50 60];    % number of voxels in ROI1 and ROI2
mVoxs = min(nVoxs);
nRuns = 2;          % number of runs. MVPD requires at least 2 runs.
sigma = 1;          % std of Noise in ROI1 and ROI2

% Optional parameters
opt.method = 'svd_exvar'; % reduce dimensions based on percent variance
opt.percentage = 90; % percent variance cut-off
opt.meancorrection = 1; % for dimension reduction (overwritten for pca_fc in data2mvpd_lprd_fc.m)
opt.regularisation = 10.^(-1:0.2:3); % regularisation parameter for ridge regression.
opt.zscore = 0; % Whether to Z-score timeseries for dCor (need to turn off for Examples 5-6)
opt.nRandomisation = 20; % Needs to be >1 to get more stable estimates of null distributions (though will trigger parfor if >1)

rng('default')

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
        Y{s}{r} = X{s}{r}*T;
        X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
        Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
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
        Y{s}{r} = X{s}{r}*T;
        X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
        Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
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
cc = 0;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Y{s}{r} = X{s}{r}*T;
        X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
        Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
    end
end
vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example3.png'),'png')

%% Fourth example: run-dependent linear mapping (within-run measures work)
fnam = 'Run-dependent linear connectivity';
cc = 0;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
% C = cov(rand(nVoxs(1))); %we can use a general cov matrix here
X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns
        X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        T = randn(nVoxs);% changed this from rand(nVoxs)-.5
        Y{s}{r} = X{s}{r}*T;
        X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
        Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
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
cc = 0.5; % cannot be zero now
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
%C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
X = {}; Y = {};
%sigma = 0.5; % Need to reduce noise?
%opt.zscore = 0; % Need to turn off z-scoring of dCor (or turn off for all examples above? - results same...)
for s=1:nSubj
    for r=1:nRuns
       X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       Y{s}{r} = abs(X{s}{r}*T);
%        Y{s}{r} = (X{s}{r}*T).^2;
        X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
        Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
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
        Y{s}{r} = X{s}{r}*T;
        X{s}{r} = X{s}{r} + sigma*randn(nTime,nVoxs(1)); % Add independent measurement noise
        Y{s}{r} = Y{s}{r} + sigma*randn(nTime,nVoxs(2));
        Y{s}{r} = Y{s}{r} + 5*sigma*repmat(randn(nTime,1),1,nVoxs(2));
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
        Y{s}{r} = X{s}{r}*T + sigma*randn(nTime,nVoxs(2));
%        X{g}{r} = X{g}{r} + sigma*randn(nTime,nVoxs(1));  % independent noise
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

meanvl(1) = mean(MVconn.RCA - MVconn_null.RCA);
spread(1) = std(MVconn.RCA - MVconn_null.RCA);
meanvl(2) = mean(sMVconn.RCA - sMVconn_null.RCA);
spread(2) = std(sMVconn.RCA - sMVconn_null.RCA);
figure('Color','w')
c = categorical({'1 RCA-orig','2 RCA-pooled'});
bar(c,meanvl,.5,'FaceColor',[0.75,0.75,0.75])
hold on
errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
ylabel('Mean +/- SD')
saveas(gcf,fullfile('Graphics','mvcon_example7.png'),'png')

%% Eighth example: multivariate lagged interaction (where lagged MV-Conn metric works better)
fnam = 'Multivariate lagged interaction';

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

% functional mapping (from ROI1 to ROI2)                 
T=randn(nVoxs);

cc = 0;
Ca = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI

X = {}; Y = {};
for s=1:nSubj
    for r=1:nRuns

        % lagged multivariate interaction
        X{s}{r}=mvnrnd(zeros(nTime+delaYin,nVoxs(1)),Ca);
        Y{s}{r}=(X{s}{r}(delaYin+1:end,:)*T)+sigma*randn(nTime,nVoxs(2));
        X{s}{r}=X{s}{r}(1:nTime,:)+sigma*randn(nTime,nVoxs(1));

    end
end

vis = [1 2 1 2];
[MVconn,MVconn_null] = computeMVconn(X,Y,opt);
plotmv(fnam,T,Ca,X,Y,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example8.png'),'png')

%% Nineth example: model-based RCA gives different results to direct RCA
clear;clc
fnam = 'model-based vs. direct RCA';
% the idea here is that the correlation of RDMs to one model can be very
% different from the correlation of RDMs.
k = 12;
ndissim = nchoosek(k,2);
T = @(theta) kron([cosd(theta),sind(theta);-sind(theta),cosd(theta)],eye(ndissim/2)); % rotation in the dissimilarity space
model = zscore(pdist(randn(k)))';
ang = 70;% this needs to be > 45 for the example to work
ds_roi1 = T(ang)*model;
ds_roi2 = T(-ang)*model;
rc_model_roi1 = corr(ds_roi1,model);
rc_model_roi2 = corr(ds_roi2,model);
rc_direct = corr(ds_roi1,ds_roi2);
c = categorical({'1 RCA-model2ROI1','2 RCA-model2ROI2','3 RCA-ROI12ROI2'});
figure('Color','w')
bar(c,[rc_model_roi1,rc_model_roi2,rc_direct],.5,'FaceColor',[0.75,0.75,0.75])
ylim([-1 1]);ylabel('RCA')
saveas(gcf,fullfile('Graphics','mvcon_example9.png'),'png')

return
