
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

CV = 0; % Whether cross-validate or not (produces different measures in case of fMRI)
win = {[1:121] [141:261]}; % indices of windows for fMRI (if use CV)
Nscram = 20;

cd(bwd)

wd = fullfile(bwd,'Data','fMRI');
try mkdir(wd); end
cd(wd)

fnames = dir('CC*');
Nsubj = length(fnames)

% Code written flexibly below so can parfor across longest dimension (subjects, ROIs or scramblings)
NumWorkers = Nsubj;
%NumWorkers = Nroi1;
%NumWorkers = Nscram+1;
P = cbupool(NumWorkers);
P.SubmitArguments = sprintf('--ntasks=%d --mem-per-cpu=4G --time=72:00:00',NumWorkers);
parpool(P,NumWorkers)

minPC = cell(Nsubj,1); allFC = minPC; wc_bc = minPC;
tic
%for subj = 1:Nsubj
parfor subj = 1:Nsubj
    fprintf('%d %s\n',subj,fnames(subj).name)
    
    minPC{subj} = cell(1,Nroi1); allFC{subj} = minPC{subj}; wc_bc{subj} = minPC{subj};
    % load data
    tmp = rikload(fnames(subj).name);
    
    % Z-score data (needed if generating phase-scrambled)
    nv = NaN; dat = {};
    for roi = rois
        dat{roi} = zscore(tmp.dat{roi});
        nv(roi) = size(dat{roi},2);
    end
    tmp = [];
    
    %parfor r1 = 1:Nroi1
    for r1 = 1:Nroi1
        fprintf('ROI %d\n',r1)
        
        % Get data for ROI1
        X = {};
        if ~CV
            X = {dat{rois1(r1)}};
        else
            for run = 1:length(win)
                X{run} = dat{rois1(r1)}(win{run},:);
            end
        end
     
        opt2 = opt; opt2.percentage = 95; opt2.meancorrection = 0;
        tmp = dimreduction(X{1},'svd_exvar',opt2);
        minPC{subj}{r1} = size(tmp,2); % Just take first run for this
        
        FC = cell(1,Nscram+1);   
        fprintf('Iteration: ')
%       parfor ps = 1:(Nscram+1) % ps==1 corresponds to true data
        for ps = 1:(Nscram+1) % ps==1 corresponds to true data
            
            fprintf('%d.',ps-1)

            fc = nan(1,Nroi2); fc_svd = fc; cca = fc;
            mvpd = fc; lprd = fc; dcor = fc; rc_w = fc; rc_b = fc; 
             
            for r2 = 1:Nroi2               
                % get data for ROI2               
                if ~CV
                    Y = {dat{rois2(r2)}};
                    if ps > 1
                        Y{1} = phase_rand(Y{1}, 1, 1); % Third argument means common scrambling across voxels
                    end
                    
                    [fc(r2),fc_svd(r2),cca(r2)] = data2UVFC(X,Y);
                else
                    Y = {};
                    for run = 1:length(win)
                        Y{run} = dat{rois2(r2)}(win{run},:);
                        if ps > 1
                            Y{run} = phase_rand(Y{run}, 1, 1); % Third argument means common scrambling across voxels
                        end
                    end
                    
                    [mvpd(r2),lprd(r2),fc(r2),fc_svd(r2),cca(r2)] = data2mvpd_lprd_fc(X,Y,opt);
                end
                [dcor(r2),~] = data2dCor(X,Y);
                [rc_w(r2),rc_b(r2)] = data2rc(X,Y,'Correlation');
            end
            
            if ~CV
                FC{ps} = struct('fc',fc,'fc_svd',fc_svd,'dcor',dcor,'rc_w',rc_w,'cca',cca);
            else
                FC{ps} = struct('fc',fc,'fc_svd',fc_svd,'mvpd',mvpd,'lprd',lprd,'dcor',dcor,'rc_w',rc_w,'rc_b',rc_b,'cca',cca);
            end
        end
        fprintf('\n')
        
        fs = fieldnames(FC{1});
        pFC = nan(length(fs),Nscram+1,Nroi2); nFC = nan(length(fs),Nroi2);
        for f=1:length(fs)
            %figure,imagesc(FC{1}.(fs{f}))
            for ps = 1:(Nscram+1)
                pFC(f,ps,:) = FC{ps}.(fs{f});
            end
            if Nscram > 1 % Care: normalisation below requires reaosnably big Nscram!
                for r2 = 1:Nroi2
                    nFC(f,r2) = (pFC(f,1,r2) - mean(pFC(f,2:end,r2))) / std(pFC(f,2:end,r2));
                end
            else
                nFC(f,:) = squeeze(pFC(f,1,:));
            end
            %figure,imagesc(nFC)
            
            % Calculate homologous minus mean on non-homologous connectivity (assuming ROIs organised that way)
            wc_bc{subj}{r1}(f) = nFC(f,r1) - mean(nFC(f,setdiff([1:Nroi2],r1)));
        end
        
        allFC{subj}{r1} = FC;       
    end
end
toc

res = []; nPC = [];
for subj = 1:Nsubj
    for r1 = 1:Nroi1
        nPC(subj,r1) = minPC{subj}{r1};        
        res(subj,r1,:) = wc_bc{subj}{r1};
    end
end

save(sprintf('fMRI_Results_CV%d',CV),'res','allFC','Nscram','nPC')

cd(bwd)

return

% Mean connections
meanROI= squeeze(mean(res,2));
fs = fieldnames(allFC{1}{1}{1});
for f=1:length(fs)
    fprintf('%s\t',fs{f})
end
fprintf('\n')
for f=1:length(fs)
    fprintf('%3.2f\t',mean(meanROI(:,f))./std(meanROI(:,f)))
    %mean(meanROI)
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



% Mean across subjects, count best method per connection 
cb = [];
meanSubj = squeeze(mean(res,1)); 
if Nscram == 0
    meanSubj = meanSubj./squeeze(std(res)); %normalise across subjects rather than phase-scramblings
end
cn = zeros(1,length(fs));
for r=1:size(meanSubj,1)
    [~,f] = max(meanSubj(r,:));
    cn(f) = cn(f)+1;
    cb(r) = f;
end
for f=1:length(fs)
    fprintf('%s\t',fs{f})
end
fprintf('\n')
for f=1:length(fs)
    fprintf('%d\t',cn(f))
end
fprintf('\n')


% Are connections where MV does better ones with more PCs:
mPC = zeros(1,length(fs));
for f=1:length(fs)
    Best = find(cb==f);
    if ~isempty(Best)
        mPC(f) = mean(mean(nPC(:,Best)));
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


