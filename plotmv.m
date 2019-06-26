function [fc, uvpd,mvpd,dcor_u,dcor, rca] = plotmv(fnam,T,Ya,Yb,opt,vis);

figure('name',fnam,'Color','w');
% plot the transformation
subplot(2,2,1)
imagesc(T)
%colormap jet
colormap gray
caxis([-1 1])
colorbar
title('Functional mapping')
xlabel('Voxels in ROI2')
ylabel('Voxels in ROI1')
set(gca,'XTick',[]); set(gca,'YTick',[])

if (vis==1)
    
    %plot the time series for the ROI1 and run1 
    v1 = 1; v2 = round(size(Ya{1}{1},2)/2)+1;
    maxY = max([Ya{1}{1}(:,v1); Ya{1}{1}(:,v2); Yb{1}{1}(:,v1); Yb{1}{1}(:,v2)]);
    minY = min([Ya{1}{1}(:,v1); Ya{1}{1}(:,v2); Yb{1}{1}(:,v1); Yb{1}{1}(:,v2)]);

    C=floor(1000*corrcoef(Ya{1}{1}(:,v1),Ya{1}{1}(:,v2)))/1000;
    subplot(2,2,3)
    plot(Ya{1}{1}(1:50,v1),':k','LineWidth',2)
    hold on
    plot(Ya{1}{1}(1:50,v2),'k','LineWidth',2)
    axis([1 50 minY maxY])
    title(strcat('ROI1, corr=',num2str(C(1,2))))
    xlabel('Time')
    ylabel('a.u.')
    legend('voxel1','voxel2','Location','best')

    % plot the time series for the ROI2 and run1
    C=floor(1000*corrcoef(Yb{1}{1}(:,1),Yb{1}{1}(:,v2)))/1000;
    subplot(2,2,4)
    plot(Yb{1}{1}(1:50,v1),':k','LineWidth',2)
    hold on
    plot(Yb{1}{1}(1:50,v2),'k','LineWidth',2)
    axis([1 50 minY maxY])
    title(strcat('ROI2, corr=',num2str(C(1,2))))
    xlabel('Time')
    ylabel('a.u.')
    legend('voxel1','voxel2','Location','best')
    
end

for g=1:length(Ya)
    [mvpd(g),uvpd(g),fc(g)] = data2mvpd(Ya{g},Yb{g},opt); 
    [dcor(g),dcor_u(g)] = data2dCor(Ya{g},Yb{g});
    [rc(g),~] = data2rc(Ya{g},Yb{g},'correlation');
end

% plot the performance
subplot(2,2,2)
c = categorical({'1 Pearson','2 UVPD','3 MVPD','4 UVdCor','5 dCor','6 RCA'});
bar(c,[mean(fc),mean(uvpd),mean(mvpd),mean(dcor_u),mean(dcor),mean(rc)],'FaceColor',[0.5,0.5,0.5])
hold
errorbar(c,[mean(fc),mean(uvpd),mean(mvpd),mean(dcor_u),mean(dcor),mean(rc)],[std(fc),std(uvpd),std(mvpd),std(dcor_u),std(dcor),std(rc)],'ko','MarkerSize',1,'CapSize',15)
temp = get(gca,'YLim');set(gca,'YLim',[temp(1)-.1,temp(2)+.1])
title('Performance')

