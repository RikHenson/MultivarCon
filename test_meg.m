

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
opt.zscore = 0; % Whether to Z-score timeseries for dCor
opt.segleng=100;    % length (in bins) of the segment used for the cross-spectrum computation,
%opt.freqbins=[30:70]+1; % low gamma band, i.e. from 30 Hz to 70 Hz
%opt.freqbins=[8:12]+1; % alpha (+1 because of way code works)
%opt.freqbins=[4:8]+1; % theta (+1 because of way code works)
%opt.freqbins=[6:48]+1; % theta,alpha,beta,low gamma
opt.freqbins=[13:30]+1; % beta

% rois1 = [7]; % L Precentral (motor) 
% rois2 = rois1+48; rois2 = [rois2 setdiff([1:48]+48,rois2)]; % All R Hemi
rois1 = [1:48];     % All HOA cortical ROIs on Left
rois2 = rois1+48;   % All HOA cortical ROIs on Right
%rois1 = [7 45]; % L Precentral (motor) and Heschl (auditory)
%rois2 = [55 93]; % R Precentral (motor) and Heschl (auditory) 
rois = unique([rois1 rois2]);
Nroi1 = length(rois1)
Nroi2 = length(rois2)

Nscram = 20;
Orthog = 0; % whether to orthogonalise (remove 0-lag) (only applies to MEG below actually)

Nsubj = 20;
NumWorkers = Nsubj;
%NumWorkers = Nroi1;
%NumWorkers = max(Nscram);
P = cbupool(NumWorkers);
P.SubmitArguments = sprintf('--ntasks=%d --mem-per-cpu=4G --time=72:00:00',NumWorkers);
parpool(P,NumWorkers)

cd(bwd)

wd = fullfile(bwd,'Data','MEG');
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
    %        parfor ps = 1:(Nscram+1) % ps==1 corresponds to true data
    minPC{subj} = nan(Nroi1,1);
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
        %            parfor r1 = 1:Nroi1
        for r1 = 1:Nroi1
            fprintf('%d.',r1)
            
            % Get data for ROI1
            X = {dat{rois1(r1)}};
            opt2 = opt; opt2.percentage = 95; opt2.meancorrection = 0;
            rX = dimreduction(X{1},'svd_exvar',opt2);
            minPC{subj}(r1) = size(rX,2); % Just take first run for this
            opt2.number = minPC{subj}(r1);
            
            if ps == 1 % won't work if ps loop is "parfor" - needs to be "for"
                X{1} = rX;
            else
                X{1} = rX(:,1:minPC{subj}(r1));
            end
             
            for r2 = 1:length(rois2)
%                for r2 = [r1 rem(r1,Nroi1)+1];
               
                % get data for ROI2
                Y = {dat{rois2(r2)}};
                
                Y{1} = dimreduction(Y{1},'svd_ndir',opt2); % Match dimensions with r1, else cannot compare with other r2
                if Orthog
                    oX   = orthog(X{1},Y{1});
                    Y{1} = orthog(Y{1},X{1});
                    X{1} = oX;
                end
                [fc{r1}(r2),fc_svd{r1}(r2)] = data2UVFC(X,Y); % Note these 0-lag correlations not very meaningful and not comparable to lagged ones below
                [mim{r1}(r2),imcoh_svd{r1}(r2),mvlagcoh{r1}(r2),lagcoh_svd{r1}(r2)] = data2lagconn(X,Y,opt);               
            end
        end
        
        % Collect up measures into structure (in case parfor)
        FC{ps} = struct('mim',cat(1,mim{:}),'imcoh_svd',cat(1,imcoh_svd{:}),'mvlagcoh',cat(1,mvlagcoh{:}),'lagcoh_svd',cat(1,lagcoh_svd{:}),'fc',cat(1,fc{:}),'fc_svd',cat(1,fc_svd{:}));
        
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
    
save(sprintf('MEG_Results_Freq%d-%d_Nscram%d',opt.freqbins(1)-1,opt.freqbins(end)-1,Nscram),'res','allFC','Nscram','minPC')

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
    fprintf('%3.2f\t',mean(meanROI(:,f))./std(meanROI(:,f)))
    %mean(meanROI)
end
fprintf('\n')


FIG=figure('name','MEG_data','Color','w','Position',[1 1 2*560 1.5*480]); hold on
c = categorical({'ImCoh','LagCoh','MIM','MVLagCoh'});
reord = [2 4 1 3];
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
saveas(gcf,fullfile('Graphics','meg_example_Subjects.png'),'png')


% Mean across subjects, count best method per connection (only works if normalised by scrambled data)
cb = [];
if Nscram > 0
    meanSubj = squeeze(mean(res(:,:,1:4),1)); % exclude non-lagged FC measures
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

FIG=figure('name','MEG_data','Color','w','Position',[1 1 2*560 1.5*480]); hold on
c = categorical({'ImCoh','LagCoh','MIM','MVLagCoh'});
reord = [2 4 1 3];
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
saveas(gcf,fullfile('Graphics','meg_example_ROIs.png'),'png')


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

