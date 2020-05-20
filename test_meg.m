

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
%opt.freqbins=[4:48]+1; % theta,alpha,beta,low gamma
opt.freqbins=[13:30]+1; % beta (eg for just motor ROIs)

rois1 = [1:48];     % All HOA cortical ROIs on Left
rois2 = rois1+48;   % All HOA cortical ROIs on Right
%rois1 = [7 45]; % L Precentral (motor) and Heschl (auditory)
%rois2 = [55 93]; % R Precentral (motor) and Heschl (auditory) 
rois = unique([rois1 rois2]);
Nroi1 = length(rois1)
Nroi2 = length(rois2)

Nscram = 20;

cd(bwd)

wd = fullfile(bwd,'Data','MEG');
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
        X = {dat{rois1(r1)}};
        opt2 = opt; opt2.percentage = 95; opt2.meancorrection = 0;
        X{1} = dimreduction(X{1},'svd_exvar',opt2);
        minPC{subj}{r1} = size(X{1},2); % Just take first run for this
        opt2.number = minPC{subj}{r1};
        
        FC = cell(1,Nscram+1);   
        fprintf('Iteration: ')
%       parfor ps = 1:(Nscram+1) % ps==1 corresponds to true data
        for ps = 1:(Nscram+1) % ps==1 corresponds to true data
            
            fprintf('%d.',ps-1)
            
            fc = nan(1,Nroi2); fc_svd = fc;
            mvpd = fc; lprd = fc; dcor = fc; rc_w = fc; rc_b = fc; cca = fc;
            mim = fc; imcoh_svd = fc; mvlagcoh = fc; lagcoh_svd = fc;
            
            for r2 = 1:Nroi2
                
                % get data for ROI2
                Y = {dat{rois2(r2)}};
                Y{1} = dimreduction(Y{1},'svd_ndir',opt2); % Match dimensions with r1, else cannot compare with other r2
                
                if ps > 1
                    Y{1} = phase_rand(Y{1}, 1, 1); % Third argument means common scrambling across voxels
                end
                
                [fc(r2),fc_svd(r2)] = data2UVFC(X,Y); % Note these 0-lag correlations not very meaningful and not comparable to lagged ones below
                [mim(r2),imcoh_svd(r2),mvlagcoh(r2),lagcoh_svd(r2)] = data2lagconn(X,Y,opt);
            end
            
            % Collect up measures into structure (in case parfor)
            FC{ps} = struct('mim',mim,'imcoh_svd',imcoh_svd,'mvlagcoh',mvlagcoh,'lagcoh_svd',lagcoh_svd,'fc',fc,'fc_svd',fc_svd);
        end
        
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

save(sprintf('MEG_Results_Freq%d-%d_Nscram%d',opt.freqbins(1)-1,opt.freqbins(end)-1,Nscram),'res','allFC','Nscram','nPC')

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
saveas(gcf,fullfile('Graphics','meg_example_Subjects_Nscram0.png'),'png')


% Mean across subjects, count best method per connection (only works if normalised by scrambled data)
cb = [];
meanSubj = squeeze(mean(res(:,:,1:4),1)); % exclude non-lagged FC measures
if Nscram == 0
    meanSubj = meanSubj./squeeze(std(res(:,:,1:4))); %normalise across subjects rather than phase-scramblings
end
cn = zeros(1,size(meanSubj,2));
for r=1:size(meanSubj,1)
    [~,f] = max(meanSubj(r,:));
    cn(f) = cn(f)+1;
    cb(r) = f;
end
for f=1:4
    fprintf('%d\t',cn(f))
end
fprintf('\n')


% Are connections where MV does better ones with more PCs:
mPC = zeros(1,4);
for f=1:length(fs)
    Best = find(cb==f);
    if ~isempty(Best)
        mPC(f) = mean(nPC(Best));
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
saveas(gcf,fullfile('Graphics','meg_example_ROIs_Nscram0.png'),'png')


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

