
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

wc_bc = cell(Nsubj,1); allFC = wc_bc; minPC = wc_bc;
tic
%for subj = 1:Nsubj
parfor subj = 1:Nsubj
    fprintf('\n%d %s\n',subj,fnames(subj).name)
    
    % load data
    tmp = rikload(fnames(subj).name);
    
    % Z-score data (needed if generating phase-scrambled)
    nv = NaN; dat = {};
    for roi = rois
        dat{roi} = zscore(tmp.dat{roi});
        nv(roi) = size(dat{roi},2);
    end
    tmp = [];
    
    odat = dat; FC = {};
    minPC{subj} = nan(Nroi1,1);
    %     parfor ps = 1:(Nscram+1) % ps==1 corresponds to true data
    for ps = 1:(Nscram+1) % ps==1 corresponds to true data
        
        dat = {};
        fprintf('Iteration %d\n',ps-1)
        if ps > 1
            for roi = rois
                dat{roi} = phase_rand(odat{roi}, 1, 1); % Third argument means common scrambling across voxels
            end
        else
            dat = odat;
        end
        
        fc = cell(Nroi1,1); for r=1:Nroi1, fc{r} = NaN(1,Nroi1); end; fc_svd = fc;
        mvpd = fc; lprd = fc; dcor = fc; rc_w = fc; rc_b = fc; cca = fc;
        mim = fc; imcoh_svd = fc; mvlagcoh = fc; lagcoh_svd = fc;
        %         parfor r1 = 1:Nroi1
        for r1 = 1:Nroi1
            fprintf('%d.',r1)
            %warning off % Just in case in parfor loop and CCA warnings because more voxels than scans
            
            % Get data for ROI1
            if ~CV
                X = {dat{rois1(r1)}};
            else
                X = {};
                for run = 1:length(win)
                    X{run} = dat{rois1(r1)}(win{run},:);
                end
            end
            
            opt2 = opt; opt2.percentage = 95; opt2.meancorrection = 0;
            rX = dimreduction(X{1},'svd_exvar',opt2);
            minPC{subj}(r1) = size(rX,2); % Just take first run for this
            opt2.number = minPC{subj}(r1);
            
            % get data for ROI2
            for r2 = 1:length(rois2)
                %                for r2 = [r1 rem(r1,Nroi1)+1];
                
                if ~CV
                    Y = {dat{rois2(r2)}};
                    [fc{r1}(r2),fc_svd{r1}(r2),cca{r1}(r2)] = data2UVFC(X,Y);
                    %[~,~,fc{r1}(r2),fc_svd{r1}(r2)] = data2mvpd_lprd_fc(X,Y,opt);
                    [dcor{r1}(r2),~] = data2dCor(X,Y);
                    [rc_w{r1}(r2),rc_b{r1}(r2)] = data2rc(X,Y,'Correlation');
                    [rc_w{r1}(r2),~] = data2rc(X,Y,'Correlation');
                else
                    Y = {};
                    for run = 1:length(win)
                        Y{run} = dat{rois2(r2)}(win{run},:);
                    end
                    [mvpd{r1}(r2),lprd{r1}(r2),fc{r1}(r2),fc_svd{r1}(r2),cca{r1}(r2)] = data2mvpd_lprd_fc(X,Y,opt);
                    [dcor{r1}(r2),~] = data2dCor(X,Y);
                    [rc_w{r1}(r2),rc_b{r1}(r2)] = data2rc(X,Y,'Correlation');
                end
            end
        end
        
        % Collect up measures into structure (in case parfor)
        if ~CV
            FC{ps} = struct('fc',cat(1,fc{:}),'fc_svd',cat(1,fc_svd{:}),'dcor',cat(1,dcor{:}),'rc_w',cat(1,rc_w{:}),'cca',cat(1,cca{:}));
        else
            FC{ps} = struct('fc',cat(1,fc{:}),'fc_svd',cat(1,fc_svd{:}),'mvpd',cat(1,mvpd{:}),'lprd',cat(1,lprd{:}),'dcor',cat(1,dcor{:}),'rc_w',cat(1,rc_w{:}),'rc_b',cat(1,rc_b{:}),'cca',cat(1,cca{:}));
        end
        fprintf('\n')
    end
    fprintf('\n')
    
    fprintf('Min/Max Dim: %d %d\n',min(minPC{subj}), max(minPC{subj}));
    
    fs = fieldnames(FC{1});
    %         figure,
    %         for f=1:length(fs)
    %             subplot(3,2,f)
    %             imagesc(FC{1}.(fs{f}))
    %             axis square
    %             title(fs{f})
    %         end
    
    % Calculate homologous minus mean on non-homologous connectivity (assuming ROIs organised that way)
    pFC = []; nFC = [];
    for f=1:length(fs)
        %figure,imagesc(FC{1}.(fs{f}))
        for ps = 1:(Nscram+1)
            pFC(ps,:,:) = FC{ps}.(fs{f});
        end
        if Nscram > 1 % Care: normalisation below requires reaosnably big Nscram!
            for r1 = 1:Nroi1
                for r2 = 1:Nroi2
                    %for r2 = [r1 rem(r1,Nroi1)+1]
                    nFC(r1,r2) = (pFC(1,r1,r2) - mean(pFC(2:end,r1,r2))) / std(pFC(2:end,r1,r2));
                end
            end
        else
            nFC = squeeze(pFC);
        end
        %figure,imagesc(nFC)
        
        for r1 = 1:Nroi1
            if Nroi1 > 1 | Nscram > 1
                wc_bc{subj}(r1,f) = nFC(r1,r1) - mean(nFC(r1,setdiff([1:Nroi2],r1)));
                %                wc_bc{subj}(r1,f) = nFC(r1,r1) - nFC(r1,rem(r1,Nroi1)+1);
            else
                wc_bc{subj}(r1,f) = nFC(r1,r1) - mean(nFC(setdiff([1:Nroi2],r1),r1));
            end
        end
    end
    
    allFC{subj} = FC;
    mean(wc_bc{subj})
    mean(wc_bc{subj})./std(wc_bc{subj})
end
toc

res = [];
for s = 1:Nsubj
    res(s,:,:) = wc_bc{s};
end

save(sprintf('fMRI_Results_CV%d_Nscram%d',CV,Nscram),'res','allFC','Nscram','minPC')

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
title('Homology Effect (across subjects)')
saveas(gcf,fullfile('Graphics',sprintf('fmri_example_Subjects_CV%d.png',CV)),'png')


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
PCs = mean(cat(2,minPC{:})');
mPC = zeros(1,length(fs));
for f=1:length(fs)
    Best = find(cb==f);
    if ~isempty(Best)
        mPC(f) = mean(PCs(Best));
    end
end
mPC
%[mean(PCs); cb]

FIG=figure('name','fMRI_data','Color','w','Position',[1 1 2*560 1.5*480]); hold on
if CV
    c = categorical({'Pearson','Pearson-SVD','Pearson-CCA','MVPD', 'dCor','Pearson-RCA','LPRD'});
    reord = [1 2 8 3 4 5 6];
else
    c = categorical({'Pearson','Pearson-SVD','Pearson-CCA', 'dCor','Pearson-RCA'});
    reord = [1 2 5 3 4];
end
meanvl = cn(reord);
bar([1:length(reord)],meanvl,'FaceColor',[0.75,0.75,0.75])
set(gca,'XTick',[1:length(reord)],'XTickLabel',c)
ylabel('Number Best')
yyaxis right
meanvl = mPC(reord);
bar([1:length(reord)],meanvl,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
ylabel('Mean Dim')
temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
title('Homology Effect (across ROIs)')
saveas(gcf,fullfile('Graphics',sprintf('fmri_example_ROIs_CV%d.png',CV)),'png')

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


