
% Could highpass filter fMRI data

addpath('/home/rh01/GitHub/MultivarCon')

clear;clc;close all

bwd = pwd;
addpath(bwd)

%% control options
opt.fc = 1;opt.mvpd = 1;opt.dcor = 1;opt.rca = 1;
opt.method = 'svd_exvar'; % reduce dimensions based on percent variance
opt.percentage = 90; % percent variance cut-off (to match demo.m)
opt.meancorrection = 1; % for dimension reduction (overwritten for pca_fc in data2mvpd_lprd_fc.m)
opt.regularisation = 10.^(-1:0.2:3); % regularisation parameter for ridge regression.
opt.zscore = 0; % Whether to Z-score timeseries for dCor (done anyway below)

rois1 = [1:48];     % All HOA cortical ROIs on Left
rois2 = rois1+48;   % All HOA cortical ROIs on Right
rois = unique([rois1 rois2]);
Nroi1 = length(rois1)
Nroi2 = length(rois2)

CV = 1; % Whether cross-validate or not (produces different measures in case of fMRI)
win = {[1:121] [141:261]}; % indices of windows for fMRI (if use CV)
Nscram = 20;

Nsubj = 20;
NumWorkers = Nsubj;
%NumWorkers = Nroi1;
%NumWorkers = max(Nscram);
P = cbupool(NumWorkers);
P.SubmitArguments = sprintf('--ntasks=%d --mem-per-cpu=4G --time=72:00:00',NumWorkers);
parpool(P,NumWorkers)

cd(bwd)

wd = fullfile(bwd,'Data','fMRI');
try mkdir(wd); end
cd(wd)

fnames = dir('CC*');


cd(bwd)

return

% Mean connections
meanROI= squeeze(mean(res,2));
fs = fieldnames(allFC{1}{1});
for f=1:length(fs)
    fprintf('%s\t',fs{f})
end
fprintf('\n')
for f=1:length(fs)
    fprintf('M=%3.2f\t',mean(meanROI(:,f)))
end
fprintf('\n')
for f=1:length(fs)
    fprintf('T=%3.2f\t',mean(meanROI(:,f))/(std(meanROI(:,f))/sqrt(size(meanROI,1))))
end
fprintf('\n')

FIG=figure('name','fMRI_data','Color','w','Position',[1 1 2*560 1.5*480]); hold on
if CV
    c = categorical({'Pearson','Pearson-SVD','Pearson-CCA','MVPD', 'dCor','Pearson-RCA','LPRD'});
    reord = [1 2 8 3 4 5 6];
else
   c = categorical({'Pearson','Pearson-SVD','Pearson-CCA', 'dCor','Pearson-RCA'});
     reord = [1 2 5 3 4];
end
meanvl = mean(meanROI(:,reord));
spread = std(meanROI(:,reord));
bar([1:length(reord)],meanvl,'FaceColor',[0.75,0.75,0.75])
set(gca,'XTick',[1:length(reord)],'XTickLabel',c)
errorbar([1:length(reord)],meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
ylabel('Mean +/- SD')
yyaxis right
bar([1:length(reord)],meanvl./spread,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
ylabel('Mean/SD')
temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
title('Homology Effect')
saveas(gcf,fullfile('Data','fMRI',sprintf('fmri_example_CV%d.png',CV)),'png')


% Mean across subjects, count best method per connection (only works if normalised by scrambled data)
cb = [];
if Nscram > 0
    meanSubj = squeeze(mean(res,1));
    cn = zeros(1,size(meanSubj,2));
    for r=1:size(meanSubj,1)
        [~,f] = max(meanSubj(r,:));
        cn(f) = cn(f)+1;
        cb(r) = f;
    end
    for f=1:length(fs)
        fprintf('%d\t',cn(f))
    end
    fprintf('\n')
end

% Are connections where MV does better ones with more PCs:
PCs = cat(2,minPC{:})';
mPCs = mean(PCs);
wm = unique(cb)
mPC = zeros(1,7);
for w=1:length(wm)
    mPC(w) = mean(mPCs(find(cb==wm(w))));
end
mPC
%[mean(PCs); cb]

% % Binarize, ie number of connections > 0
% bres = [];
% for s = 1:size(res,1)
%     for f = 1:size(res,3)
%         bres(s,f) = length(find(res(s,:,f)>0));
%     end
% end
% bres = mean(bres);
% for f=1:length(fs)
%     fprintf('%3.1f\t',bres(f))
% end
% fprintf('\n')

return


