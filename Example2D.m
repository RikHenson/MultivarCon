cd /home/rh01/GitHub/MultivarCon
options.number = 1;
options.meancorrection = 0;

%% 1 When Pearson of Mean works

X = [
    -3  -2
    -2  -3
     3   2
     2   3
    ];

Y = X; 
X1=dimreduction(X,'svd_ndir',options);
Y1=dimreduction(Y,'svd_ndir',options);
corr(mean(X,2),mean(Y,2))
corr(X1,Y1) % X1,Y1 projections along y=x

figure,
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
xlabel('Voxel 1')
ylabel('Voxel 2')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]),
xlabel('Activity of voxel 1')
ylabel('Activity of voxel 2')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
print -dpng 'UVmean.png'


%% When Pearson of SVD works
    
Y = [-X(:,1) X(:,2)];

corr(mean(X,2),mean(Y,2))
X1=dimreduction(X,'svd_ndir',options);
Y1=dimreduction(Y,'svd_ndir',options);
corr(X1,Y1)
%corr(diff(X,1,2),sum(Y,2)) % equivalent to SVD, ie projecting along y=-x

figure
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
xlabel('Activity of voxel 1')
ylabel('Activity of voxel 2')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]),
xlabel('Activity of voxel 1')
ylabel('Activity of voxel 2')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
print -dpng 'UVsvd.png'


%% When MV needed, because data not dominated by one direction 
    
X = [
     0 -3
    -3  0
     1  1
     2  2
    ];

Y = [
     -3  -3
      2   2
     -1  -1
      1   1
];

corr(mean(X,2),mean(Y,2))
X1=dimreduction(X,'svd_ndir',options);
Y1=dimreduction(Y,'svd_ndir',options);
corr(X1,Y1) 

DX = pdist(X,'euclidean');
DY = pdist(Y,'euclidean');
corr(DX',DY')

T=(Y'*X*pinv(X'*X));
Y2=(T*X')';

figure('position',[300 300 1200 400]);
subplot(1,3,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
xlabel('Activity of voxel 1')
ylabel('Activity of voxel 2')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
subplot(1,3,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]),
xlabel('Activity of voxel 1')
ylabel('Activity of voxel 2')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
subplot(1,3,3); title('Transformed ROI1/estimated ROI2')
hold on, for t=1:size(X,1); t = text(Y2(t,1),Y2(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4])%, set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]),
xlabel('First combination of voxel activities')
ylabel('Second combination of voxel activities')
line([0 0], [-4 4],'Color','red','LineStyle',':','LineWidth',0.5)
line([-4 4], [0 0],'Color','red','LineStyle',':','LineWidth',0.5)
print -dpng 'MV1.png'

% RDM representation
figure,
subplot(1,2,1); imagesc(squareform(DX)), colormap('gray'), set(gca,'XTick',[1:4]), set(gca,'YTick',[1:4]), axis square, title('ROI1')
xlabel('Voxel 1')
ylabel('Voxel 2')
subplot(1,2,2); imagesc(squareform(DY)), colormap('gray'), set(gca,'XTick',[1:4]), set(gca,'YTick',[1:4]), axis square, title('ROI2')
xlabel('Voxel 1')
ylabel('Voxel 2')
