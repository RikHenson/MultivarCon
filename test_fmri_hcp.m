% test_fmri_hcp.m
% compute the MV connectivity measures for two ROIs in the HCP dataset
% dependencies: you need the cifti toolbox 



warning off
subj     = 'SUBJID'; % this is the subject number for the HCP data
sessions = {'1_LR' '1_RL' '2_LR' '2_RL'};
roi1 = 120;roi2 = 180+roi1;
%% control options
opt.fc = 1;opt.mvpd = 1;opt.dcor = 1;opt.rca = 1;
opt.method = 'svd_exvar'; % reduce dimensions based on percent variance
opt.percentage = 90; % percent variance cut-off
opt.meancorrection = 1;
opt.regularisation = 10.^(-1:0.2:3); % regularisation parameter for ridge regression.
opt.zscore = 0; % Whether to Z-score timeseries for dCor
opt.nRandomisation = 20;
%% load atlas parcellations
indices = defineROIs_Glasser(subj);
%% load the rest data for the full brain
sessions = {'1_LR' '1_RL' '2_LR' '2_RL'};
X = cell(numel(sessions),1);Y = X;
for s = 1:length(sessions)
    session = sessions{s};
    disp(session)
    datadir = ['/vols/Scratch/HCP/rfMRI/subjects1200/' subj '/MNINonLinear/Results/rfMRI_REST' session];
    fname   = [datadir '/rfMRI_REST' session '_Atlas_MSMAll_hp2000_clean.dtseries.nii' ];
    
    [cifti, xml] = cifti_open(fname);
    data = double(cifti.cdata);
    % choose data from rois 1 & 2
    X{s} = data(indices{roi1},:)';
    Y{s} = data(indices{roi2},:)';  
end

%% compute connectivity for the two regions
fprintf('computing raw connectivity estimates...\n')
[mvpd,lprd,fc,fc_svd] = data2mvpd_lprd_fc(X,Y,opt);
[~,~,fc_cca]=data2UVFC(X,Y);
dcor = data2dCor(X,Y);
rc = data2rc(X,Y,'correlation');

%% shuffle the time points (consistently for all voxels) in each roi
for sess=1:4
    pats_null_x{sess} = phase_rand(X{sess},opt.nRandomisation,1);
    pats_null_y{sess} = phase_rand(Y{sess},opt.nRandomisation,1);
end

%% compute the null values
bmvpd= {}; blprd = {}; bfc = {}; bfc_svd = {}; bdcor = {}; 
brc = {};bfc_cca = {};
X_null = cell(numel(sessions),1);Y_null = X_null;
for iter = 1:opt.nRandomisation
    fprintf('null iteration %d out of %d\n',iter,opt.nRandomisation)
    X_null = {}; Y_null = {};
    for sess = 1:4
        X_null{sess} = squeeze(pats_null_x{sess}(:,:,iter));
        Y_null{sess} = squeeze(pats_null_y{sess}(:,:,iter));
    end
    [bmvpd{iter},blprd{iter},bfc{iter},bfc_svd{iter}] = data2mvpd_lprd_fc(X_null,Y_null,opt);
    bdcor{iter} = data2dCor(X_null,Y_null);
    brc{iter} = data2rc(X_null,Y_null,'Correlation');
    [~,~,bfc_cca{iter}]=data2UVFC(X,Y);
end
toc
resdir = '/Users/hnili/Desktop/HC/MDconValues';
if ~exist(resdir),mkdir(resdir);end
save([resdir '/' subj '.mat'],'mvpd','bmvpd','lprd','blprd',...
    'dcor','bdcor','rc','brc','fc','bfc','fc_svd','bfc_svd','fc_cca',...
    'bfc_cca');

%% visualisation of the results
cd(resdir)
fnames = dir('*.mat');
N = size(fnames,1);

for s=1:N
    disp(s)
    load(fnames(s).name);
    values = [fc,fc_svd,fc_cca,mvpd,dcor,rc,lprd];
    bvalues = [cell2mat(bfc)',cell2mat(bfc_svd)',cell2mat(bfc_cca)',...
        cell2mat(bmvpd)',cell2mat(bdcor)',cell2mat(brc)',cell2mat(blprd)'];
    zs(s,:) = (values-mean(bvalues))./nanstd(bvalues);
    meanVals(s,:) = values;
end

%% display results (normalised across and within subjects)
close all;figw(1)
meth = {'Pearson','Pearson-SVD','Pearson-CCA','MVPD', 'dCor','Pearson-RCA'};
meanVals = zs(:,1:6);
meanvl = nanmean(meanVals);
spread = nanstd(meanVals);
bar(1:numel(meth),meanvl,'FaceColor',[0.75,0.75,0.75])
set(gca,'XTick',1:numel(meth),'XTickLabel',meth)
hold on
errorbar(1:numel(meth),meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
ylabel('Mean +/- SD')
yyaxis right
bar([1:numel(meth)],meanvl./spread,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
ylabel('Mean/SD')
temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
xtickangle(45)

