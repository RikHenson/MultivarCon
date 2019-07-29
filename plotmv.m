function [] = plotmv(fnam,T,C,X,Y,MVconn,MVconn_null,vis)
% visualisation function for the results of computeMVconn.m
FIG=figure('name',fnam,'Color','w','Position',[1 1 2*560 1.5*480]);

% plot the ROI covariance
subplot(3,2,1)
imagesc(C)
%colormap jet
colormap gray
axis square
caxis([-1 1])
colorbar
title('A. ROI1 voxel covariance')
xlabel('Voxels in ROI1')
ylabel('Voxels in ROI1')
set(gca,'XTick',[]); set(gca,'YTick',[])

% plot the transformation
subplot(3,2,2)
imagesc(T)
%colormap jet
colormap gray
%axis square
caxis([-1 1])
colorbar
title('B. Functional mapping')
xlabel('Voxels in ROI2')
ylabel('Voxels in ROI1')
set(gca,'XTick',[]); set(gca,'YTick',[])

nTime = min(size(X{1}{1},1),50);
%plot the timeseries for first nTime timepoints in run 1...
if length(vis) == 4  % assume passing 2 voxel indices
    
    %for just one voxel from ROI1...
    maxY = max([X{1}{1}(:,vis(1)); X{1}{1}(:,vis(2)); Y{1}{1}(:,vis(3)); Y{1}{1}(:,vis(4))]);
    minY = min([X{1}{1}(:,vis(1)); X{1}{1}(:,vis(2)); Y{1}{1}(:,vis(3)); Y{1}{1}(:,vis(4))]);
    
    C=floor(1000*corrcoef(X{1}{1}(:,vis(1)),X{1}{1}(:,vis(2))))/1000;
    
    subplot(3,2,3)
    plot(X{1}{1}(1:nTime,vis(1)),':k','LineWidth',2)
    hold on
    plot(X{1}{1}(1:nTime,vis(2)),'k','LineWidth',2)
    axis([1 nTime minY maxY])
    title(strcat('C. ROI1 data, corr=',num2str(C(1,2))))
    xlabel('Time or Trial')
    ylabel('a.u.')
    legend(sprintf('voxel%d',vis(1)),sprintf('voxel%d',vis(2)),'Location','best')
    
    % ...and one voxel from ROI 2
    C=floor(1000*corrcoef(Y{1}{1}(:,1),Y{1}{1}(:,vis(4))))/1000;
    subplot(3,2,4)
    plot(Y{1}{1}(1:nTime,vis(3)),':k','LineWidth',2)
    hold on
    plot(Y{1}{1}(1:nTime,vis(4)),'k','LineWidth',2)
    axis([1 nTime minY maxY])
    title(strcat('D. ROI2 data, corr=',num2str(C(1,2))))
    xlabel('Time or Trial')
    ylabel('a.u.')
    legend(sprintf('voxel%d',vis(3)),sprintf('voxel%d',vis(4)),'Location','best')
    
elseif vis == 2
    
    %...or all voxels from ROI1
    %    maxY = max([X{1}{1}(:); Y{1}{1}(:)]);
    %    minY = min([X{1}{1}(:); Y{1}{1}(:)]);
    maxY = max(X{1}{1}(:));
    minY = min(X{1}{1}(:));
    
    subplot(3,2,3)
    imagesc(X{1}{1}(1:nTime,:)')
    caxis([minY maxY])
    colorbar
    xlabel('Time')
    ylabel('voxel')
    title('C. ROI1 data')
    
    %...and all voxels from ROI2
    maxY = max(Y{1}{1}(:));
    minY = min(Y{1}{1}(:));
    
    subplot(3,2,4)
    imagesc(Y{1}{1}(1:nTime,:)')
    caxis([minY maxY])
    colorbar
    xlabel('Time')
    ylabel('voxel')
    title('D. ROI2 data')
end

% plot the absolute performance
if ~isfield(MVconn,'MIM')
    subplot(3,2,5), hold on
    c = categorical({'1 Pearson','2 Pearson-PCA','3 MVPD','4 GOF','5 dCor','6 RCA'});
    allval = [MVconn.FC MVconn.FCPC MVconn.MVPD MVconn.GOF MVconn.dCor MVconn.RCA];
    meanvl = mean(allval);
    spread = std(allval);
    % spread = iqr(allval);
    bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
    errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
    ylabel('Mean +/- SD')
    yyaxis right
    bar(c,meanvl./spread,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
    ylabel('Mean/SD')
    temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
    title('E. Raw Performance')
else
    close(FIG)
    FIG=figure('name',fnam,'Color','w','Position',[1 1 2*560 1*480]);
    subplot(1,2,1), hold on
    c = categorical({'1 ImCoh','2 ImCoh-PCA','3 MIM'});
    allval = [MVconn.ImCoh MVconn_null.ImCohPC MVconn.MIM];
    meanvl = mean(allval);
    spread = std(allval);
    % spread = iqr(allval);
    bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
    errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
    ylabel('Mean +/- SD')
    yyaxis right
    bar(c,meanvl./spread,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
    ylabel('Mean/SD')
    temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
    title('A. Raw Performance')
end


% Calculate connectivity when X and Y independent random noise (since
% some connectivity measures, eg Dcor, not bounded by 0 or -1)
if length(X)<20
    warning('Insufficient subjects (<20) to estimate baseline error')
else
    if ~isfield(MVconn,'MIM')
        subplot(3,2,6), hold on
        allnul = [MVconn_null.FC MVconn_null.FCPC MVconn_null.MVPD MVconn_null.GOF MVconn_null.dCor MVconn_null.RCA];
        % meanvl = meanvl - mean(allnul);
        % spread = sqrt(spread.^2 + var(allnul));
        meanvl = mean(allval - allnul);
        spread = std(allval - allnul);
        bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
        errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
        ylabel('Mean +/- SD')
        yyaxis right
        bar(c,meanvl./spread,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
        ylabel('Mean/SD')
        temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
        title('F. Normalised Performance')
    else
        subplot(1,2,2), hold on
        allnul = [MVconn_null.ImCoh MVconn_null.ImCohPC MVconn_null.MIM];
        % meanvl = meanvl - mean(allnul);
        % spread = sqrt(spread.^2 + var(allnul));
        meanvl = mean(allval - allnul);
        spread = std(allval - allnul);
        bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
        errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
        ylabel('Mean +/- SD')
        yyaxis right
        bar(c,meanvl./spread,0.4,'FaceColor',[0 0 1],'FaceAlpha',0.3)
        ylabel('Mean/SD')
        temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
        title('B. Normalised Performance')        
    end
end
