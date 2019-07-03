function [] = plotmv(fnam,T,C,Ya,Yb,MVconn,MVconn_null,vis)
% visualisation function for the results of computeMVconn.m

figure('name',fnam,'Color','w','Position',[1 1 2*560 1.5*480]);

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

nTime = min(size(Ya{1}{1},1),50);

%plot the timeseries for first nTime timepoints in run 1... 
if length(vis) == 4  % assume passing 2 voxel indices

    %...for just one voxel from ROI1
    maxY = max([Ya{1}{1}(:,vis(1)); Ya{1}{1}(:,vis(2)); Yb{1}{1}(:,vis(3)); Yb{1}{1}(:,vis(4))]);
    minY = min([Ya{1}{1}(:,vis(1)); Ya{1}{1}(:,vis(2)); Yb{1}{1}(:,vis(3)); Yb{1}{1}(:,vis(4))]);
    
    C=floor(1000*corrcoef(Ya{1}{1}(:,vis(1)),Ya{1}{1}(:,vis(2))))/1000;
    
    subplot(3,2,3)
    plot(Ya{1}{1}(1:nTime,vis(1)),':k','LineWidth',2)
    hold on
    plot(Ya{1}{1}(1:nTime,vis(2)),'k','LineWidth',2)
    axis([1 nTime minY maxY])
    title(strcat('C. ROI1 data, corr=',num2str(C(1,2))))
    xlabel('Time or Trial')
    ylabel('a.u.')
    legend(sprintf('voxel%d',vis(1)),sprintf('voxel%d',vis(2)),'Location','best')

    % ...and one voxel from ROI 2
    C=floor(1000*corrcoef(Yb{1}{1}(:,1),Yb{1}{1}(:,vis(4))))/1000;
    subplot(3,2,4)
    plot(Yb{1}{1}(1:nTime,vis(3)),':k','LineWidth',2)
    hold on
    plot(Yb{1}{1}(1:nTime,vis(4)),'k','LineWidth',2)
    axis([1 nTime minY maxY])
    title(strcat('D. ROI2 data, corr=',num2str(C(1,2))))
    xlabel('Time or Trial')
    ylabel('a.u.')
    legend(sprintf('voxel%d',vis(3)),sprintf('voxel%d',vis(4)),'Location','best')

elseif vis == 2
    
    %...or all voxels from ROI1
%    maxY = max([Ya{1}{1}(:); Yb{1}{1}(:)]);
%    minY = min([Ya{1}{1}(:); Yb{1}{1}(:)]);
    maxY = max(Ya{1}{1}(:));
    minY = min(Ya{1}{1}(:));
     
    subplot(3,2,3)
    imagesc(Ya{1}{1}(1:nTime,:)')
    caxis([minY maxY])
    colorbar
    xlabel('Time')
    ylabel('voxel')
    title('C. ROI1 data')
 
    %...and all voxels from ROI2
    maxY = max(Yb{1}{1}(:));
    minY = min(Yb{1}{1}(:));

    subplot(3,2,4)
    imagesc(Yb{1}{1}(1:nTime,:)')
    caxis([minY maxY])
    colorbar
    xlabel('Time')
    ylabel('voxel')
    title('D. ROI2 data')
end

% plot the absolute performance
subplot(3,2,5), hold on
c = categorical({'1 Pearson','2 UVPD','3 MVPD','4 UVdCor','5 dCor','6 RCA'});
meanvl = mean([MVconn.FC MVconn.UVPD MVconn.MVPD MVconn.dCor_univar MVconn.dCor MVconn.RCA]);
spread = std([MVconn.FC MVconn.UVPD MVconn.MVPD MVconn.dCor_univar MVconn.dCor MVconn.RCA]);
% spread = iqr([MVconn.fc MVconn.uvpd MVconn.mvpd MVconn.dcor_u MVconn.dcor MVconn.rc]);
bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
title('E. Raw Performance')

% Calculate connectivity when Ya and Yb independent random noise (since
% some connectivity measures, eg Dcor, not bounded by 0 or -1)
if length(Ya)<20
    warning('Insufficient subjects (<20) to estimate baseline error')
else    
    subplot(3,2,6), hold on
    meanvl = meanvl - mean([MVconn_null.FC MVconn_null.UVPD MVconn_null.MVPD MVconn_null.dCor_univar MVconn_null.dCor MVconn_null.RCA]);
%    spread = sqrt(spread.^2 + var([MVconn_null.FC MVconn_null.UVPD MVconn_null.MVPD MVconn_null.dCor_univar MVconn_null.dCor MVconn_null.RCA]));
    spread = std([MVconn.FC MVconn.UVPD MVconn.MVPD MVconn.dCor_univar MVconn.dCor MVconn.RCA] - ...
        [MVconn_null.FC MVconn_null.UVPD MVconn_null.MVPD MVconn_null.dCor_univar MVconn_null.dCor MVconn_null.RCA]); % better, since each baseline paired with real subject?
    bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
    errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
    %bar(c,meanvl./spread,'FaceColor',[0.25,0.25,0.25])
    temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
    title('F. Normalised Performance')
end

% subplot(3,2,6), hold on
% meanvl = mean([bfc buvpd bmvpd bdcor_u bdcor brc]);
% spread = std([fc uvpd mvpd dcor_u dcor rc]);
% bar(c,meanvl,'FaceColor',[0.75,0.75,0.75])
% errorbar(c,meanvl,spread,'ko','MarkerSize',1,'CapSize',15)
% temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
% title('Normalised Performance')





