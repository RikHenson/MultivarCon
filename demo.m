close all
clear;
clc;
rng('default')

% Parameter settings
nSubj = 20;         % number of subjects (replications)
nTime = 200;        % number of time poinTimes
nVoxs = [50 60];    % number of voxels in ROI1 and ROI2
mVoxs = min(nVoxs);
nRuns = 2;          % number of runs. MVPD requires at least 2 runs.
Noise = 1;          % std of Noise in ROI2

% MVPD parameters
opt.method='pca_exvar';
opt.percentage=99;
opt.number=1;
opt.regularisation=10.^(-1:0.1:3); % regularisation parameter for ridge regression.

%% First example: positively correlated voxel activities within ROI1 (where both UniConn and MultiConn work)

fnam = 'Positively correlated voxel activities within a ROI';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nRuns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = [1 2 1 2];
%vis = 2;
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)

saveas(gcf,'Graphics/mvcon_example1.png','png')

%% Second example: anticorrelated voxel activities within ROI1 (where MulitCon work better)

fnam = 'Negatively correlated voxel activities within a ROI';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.9;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + (1-cc)*eye(nVoxs(1));
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nRuns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
saveas(gcf,'Graphics/mvcon_example2.png','png')

%% Third example: anticorrelation in ROI2 induced by the functional mapping (where MulitCon work better)

fnam = 'Negative correlations induced by the functional mapping';
%T = rand(nVoxs)-0.5;
T = zeros(nVoxs); for j=1:mVoxs/2; T(j,j)=1; T(j+mVoxs/2,j+mVoxs/2)=-1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nRuns
        Ya{g}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{g}{r} = Ya{g}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
saveas(gcf,'Graphics/mvcon_example3.png','png')

%% Fourth example: run-dependent linear mapping (within-run measures work)

fnam = 'Run-dependent linear connectivity';
cc = 0.5;
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
C = cov(rand(nVoxs(1))); %this example doesn't require a particular cov
% structure in either ROI.
Ya = {}; Yb = {};
for s=1:nSubj
    for r=1:nRuns
       Ya{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       T = randn(nVoxs);% changed this from rand(nVoxs)-.5 
       Yb{s}{r} = Ya{s}{r}*T + Noise*randn(nTime,nVoxs(2));
    end
end
vis = 2;
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example4.png'),'png')

%% Fifth example: the functional mapping is nonlinear (dCor works best)

fnam = 'Nonlinear connectivity';
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
cc = 0.5;
C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
%C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
Ya = {}; Yb = {};
Noise = 0.5;
for s=1:nSubj
    for r=1:nRuns
       Ya{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
       Yb{s}{r} = abs(Ya{s}{r}*T) + Noise*randn(nTime,nVoxs(2));
%        Yb{s}{r} = (Ya{s}{r}*T).^2 + Noise*randn(nTime,nVoxs(2));% also a good example
    end
end
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
vis = 2;
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
saveas(gcf,fullfile('Graphics','mvcon_example5.png'),'png')

%% Sixth example: the presence of structured noise in ROI2 (where RCA works best)
% do we have to zscore the data for the GOF?

fnam = 'Structured Noise in ROI2'; 
%T = rand(nVoxs);
T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end
% T = randn(nVoxs);
cc = 0.5;
%C = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc); % correlation within ROI
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + cc*eye(nVoxs(1));
% C = cov(rand(nVoxs(1)));% this can also be general. again we don't need a 
% specfic C for this
Ya = {}; Yb = {};
Noise = 1;
for s=1:nSubj
    for r=1:nRuns
        Ya{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
        Yb{s}{r} = Ya{s}{r}*T + Noise*randn(nTime,nVoxs(2));
        Yb{s}{r} = Yb{s}{r}   + 5*Noise*repmat(randn(nTime,1),1,nVoxs(2));
    end
end
vis = 2;
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
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

Ya = {}; Yb = {};
for s=1:nSubj
    for r=1:nRuns
        Ya{s}{r} = mvnrnd(zeros(nStim,nVoxs(1)),C);
        Ya{s}{r} = repmat(Ya{s}{r},nRep,1);  % Repeat same pattern 
        Yb{s}{r} = Ya{s}{r}*T + Noise*randn(nTime,nVoxs(2));
%        Ya{g}{r} = Ya{g}{r} + Noise*randn(nTime,nVoxs(1));  % independent noise
    end
end

vis = 2;
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt);
% plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
% saveas(gcf,'Graphics/mvcon_example7a.png','png')
rc_orig = mean(MVconn.RCA)-mean(MVconn_null.RCA);

sYa = {}; sYb = {};
for s=1:nSubj
    for r=1:nRuns
        for stim = 1:nStim
            sYa{s}{r}(stim,:) = mean(Ya{s}{r}(find(stimuli==stim),:));
            sYb{s}{r}(stim,:) = mean(Yb{s}{r}(find(stimuli==stim),:));
        end
    end
end

vis = 2;
%vis = [1 nVoxs(1)/2+1 1 nVoxs(2)/2+1];
[sMVconn,sMVconn_null] = computeMVconn(sYa,sYb,opt);
% plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
% saveas(gcf,'Graphics/mvcon_example7b.png','png')
rc_pooled = mean(sMVconn.RCA)-mean(sMVconn_null.RCA);
fprintf('RC from the original trial x trial RDMs: \t: %.2f\n',rc_orig)
fprintf('RC from the pooled data i.e. stim x stim RDMs: \t: %.2f\n',rc_pooled)



%% Eighth example: closed-loop (where lagged MultiCon work better)

fnam = 'Closed-loop between subpopulations';
 % functional mapping between the two regions

% Parameter settings
nTime = 15360;      % number of time points, e.g. 15360 is equivalent to 60 seconds 
                    % acquisition if the resolution is 256 Hz
nVoxs = [6 3];      % number of voxels in ROI1 and ROI2

% MIM parameters
opt.segleng=256;    % length (in bins) of the segment used for the cross-spectrum computation,
                    % e.g. 256 is equivalent to 1 second segment with a resolution of 256 Hz
opt.freqbins=31:71; % low gamma band, i.e. from 30 Hz to 70 Hz

delaybin= 10;     % delay in bins of the interaction between the regions,
                  % e.g. if the resolution is 256 Hz and the delaybin is 10, 
                  % the actual delay between a and b is delaybin/(256 Hz)= 40ms 
Noise=0.5;
                  
% correlated (mixed) noise
Mixnoise=randn(nTime,sum(nVoxs));
Mixnoisea=Mixnoise*randn(sum(nVoxs),nVoxs(1));
Mixnoiseb=Mixnoise*randn(sum(nVoxs),nVoxs(2));

Tab=randn(nVoxs);
Tba=randn(nVoxs(2:-1:1));
cc = 0.5;
Ca = ones(nVoxs(1))*cc + eye(nVoxs(1))*(1-cc);
Cb = ones(nVoxs(2))*cc + eye(nVoxs(2))*(1-cc);


Ya = {}; Yb = {};
for g=1:nSubj
    for r=1:nRuns
        
        % simulation of lagged multivariate interaction from a population
        % in ROIa to a population in ROIb
        Ya{g}{r}=mvnrnd(zeros(nTime+delaybin,nVoxs(1)),Ca);
        %Yb{g}{r}=(Ya{g}{r}(delaybin+1:end,:)*Tab)/norm(Ya{g}{r}(delaybin+1:end,:)*Tab,'fro')+Noise*Mixnoiseb/norm(Mixnoiseb,'fro');
        %Ya{g}{r}=Ya{g}{r}(1:nTime,:)/norm(Ya{g}{r}(1:nTime,:),'fro')+Noise*Mixnoisea/norm(Mixnoisea,'fro');
        Yb{g}{r}=(Ya{g}{r}(delaybin+1:end,:)*Tab)+Noise*Mixnoiseb;
        Ya{g}{r}=Ya{g}{r}(1:nTime,:)+Noise*Mixnoisea;
        Za{g}{r}=Ya{g}{r};
        Zb{g}{r}=Yb{g}{r};

        % simulation of lagged multivariate interaction from a population
        % in ROIb to a population in ROIa
        Yb{g}{r}=mvnrnd(zeros(nTime+delaybin,nVoxs(2)),Cb);
        %Ya{g}{r}=(Yb{g}{r}(delaybin+1:end,:)*Tba)/norm(Yb{g}{r}(delaybin+1:end,:)*Tba,'fro')+Noise*Mixnoisea/norm(Mixnoisea,'fro');
        %Yb{g}{r}=Yb{g}{r}(1:nTime,:)/norm(Yb{g}{r}(1:nTime,:),'fro')+Noise*Mixnoiseb/norm(Mixnoiseb,'fro');
        Ya{g}{r}=(Yb{g}{r}(delaybin+1:end,:)*Tba)+Noise*Mixnoisea;
        Yb{g}{r}=Yb{g}{r}(1:nTime,:)+Noise*Mixnoiseb;
        Za{g}{r}=[Za{g}{r}, Ya{g}{r}];
        Zb{g}{r}=[Zb{g}{r}, Yb{g}{r}];

    end
end

vis = [1 2 1 2];
%vis = 2;
[MVconn,MVconn_null] = computeMVconn(Za,Zb,opt);
plotmv(fnam,Tab,Ca,Za,Zb,MVconn,MVconn_null,vis)
saveas(gcf,'Graphics/mvcon_example8.png','png')

return

