clear;clc;close all
% Parameter settings
nSubj = 20;         % number of subjects (replications)
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
sigma = 1;
rng('default')

%% simulate
nTime = sum(nVoxs)*[1:0.5:10];
T = randn(nVoxs); cc = 0; % Demo 3
C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + (1-cc)*eye(nVoxs(1));
X = {}; Y = {};
for g=1:length(nTime)
    clear X Y
    for s=1:nSubj
        for r=1:nRuns
            X{s}{r} = mvnrnd(zeros(nTime(g),nVoxs(1)),C);
            Y{s}{r} = X{s}{r}*T;
            X{s}{r} = X{s}{r} + sigma*randn(nTime(g),nVoxs(1)); % Add independent measurement noise
            Y{s}{r} = Y{s}{r} + sigma*randn(nTime(g),nVoxs(2));
        end
    end
    [MVconn(g),MVconn_null(g)] = computeMVconn(X,Y,opt);
    disp(g)
end


%% plot the shaded error bar
methods = {'FC','FCSVD','FCCCA','MVPD','dCor','RCA','LPRD'};
for g=1:length(nTime)
    for meth = 1:numel(methods)
        allval = MVconn(g).(methods{meth});
        allnull = MVconn_null(g).(methods{meth});
        y(g,meth) = mean(allval-allnull);
        se(g,meth) = std(allval-allnull);
    end
end

figure;hold on
nTime=[1:0.5:10];
s1=shadedErrorBar(nTime,y(:,1),se(:,1),'k');
s2=shadedErrorBar(nTime,y(:,2),se(:,2),'b');
s3=shadedErrorBar(nTime,y(:,3),se(:,3),'m');
s4=shadedErrorBar(nTime,y(:,4),se(:,4),'c');
s5=shadedErrorBar(nTime,y(:,5),se(:,5),'g');
s6=shadedErrorBar(nTime,y(:,6),se(:,6),'y');
s7=shadedErrorBar(nTime,y(:,7),se(:,7),'r');
legend([s1.mainLine s2.mainLine s3.mainLine s4.mainLine s5.mainLine s6.mainLine s7.mainLine],methods,'Location','NorthWest');
xlabel('#time-points/total # voxels in the two ROIs')
ylabel('Normalised Performance')
set(gca,'FontSize',18)
axis([1 10 -0.1 1])
saveas(gcf,fullfile('Graphics','mvcon_example3_timepoints.png'),'png')

