clear;clc;close all
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
nGammas = 20;

% define function for convex combination
get_data = @(s,n,gamma) gamma*s/norm(s,'fro')+(1-gamma)*n/norm(n,'fro');
rng('default')

%% simulate
gammas = linspace(0,1,nGammas);
%T = zeros(nVoxs); for j=1:mVoxs; T(j,j)=1; end; cc = 0.9; % Demo 2
T = randn(nVoxs); cc = 0; % Demo 3

C = kron([cc -cc; -cc cc],ones(nVoxs(1)/2)) + (1-cc)*eye(nVoxs(1));
X = {}; Y = {};
for g=1:nGammas
    for s=1:nSubj
        for r=1:nRuns
            X{s}{r} = mvnrnd(zeros(nTime,nVoxs(1)),C);
            Y{s}{r} = X{s}{r}*T;
            X{s}{r} = get_data(X{s}{r},randn(nTime,nVoxs(1)),gammas(g));
            Y{s}{r} = get_data(Y{s}{r},randn(nTime,nVoxs(2)),gammas(g));
        end
    end
    [MVconn(g),MVconn_null(g)] = computeMVconn(X,Y,opt);
end


%% plot the shaded error bar
methods = {'FC','FCSVD','MVPD','LPRD','dCor','RCA'};
for g=1:nGammas
    for meth = 1:numel(methods)
        allval = MVconn(g).(methods{meth});
        allnull = MVconn_null(g).(methods{meth});
        y(g,meth) = mean(allval-allnull);
        se(g,meth) = std(allval-allnull);
    end
end
        
figure;hold on
s1=shadedErrorBar(gammas,y(:,1),se(:,1),'k');
s2=shadedErrorBar(gammas,y(:,2),se(:,2),'b');
s3=shadedErrorBar(gammas,y(:,3),se(:,3),'m');
s4=shadedErrorBar(gammas,y(:,4),se(:,4),'c');
s5=shadedErrorBar(gammas,y(:,5),se(:,5),'r');
s6=shadedErrorBar(gammas,y(:,6),se(:,6),'g');
legend([s1.mainLine s2.mainLine s3.mainLine s4.mainLine s5.mainLine s6.mainLine],methods,'Location','NorthWest');
xlabel('Signal:Noise Ratio')
ylabel('Normalised Performance')
set(gca,'FontSize',18)
set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 1 Inf])
saveas(gcf,fullfile('Graphics','mvcon_example3_noise.png'),'png')

