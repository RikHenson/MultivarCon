
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
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
print -dpng 'UVmean.png'


%% When Pearson of SVD works
    
Y = [-X(:,1) X(:,2)];

corr(mean(X,2),mean(Y,2))
X1=dimreduction(X,'svd_ndir',options);
Y1=dimreduction(Y,'svd_ndir',options);
corr(X1,Y1)
%corr(diff(X,1,2),sum(Y,2)) % equivalent to SVD, ie projecting along y=-x

figure, 
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
print -dpng 'UVsvd.png'


%% When MV needed, because data not dominated by one direction 
    
X = [
     0 -3
    -3  0
     1  1
     2  2
    ];

Y = [
     3  -3
    -3   3
    -1/2 0
     0 1/2
];

corr(mean(X,2),mean(Y,2))
X1=dimreduction(X,'svd_ndir',options);
Y1=dimreduction(Y,'svd_ndir',options);
corr(X1,Y1) 
%corr(sum(X,2),diff(Y,1,2))  % equivalent to SVD, ie projecting along y=-x

DX = pdist(X,'euclidean');
DY = pdist(Y,'euclidean');
corr(DX',DY')

figure, 
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
print -dpng 'MV1.png'

% RDM representation
figure,
subplot(1,2,1); imagesc(squareform(DX)), colormap('gray'), set(gca,'XTick',[1:4]), set(gca,'YTick',[1:4]), axis square, title('ROI1') 
subplot(1,2,2); imagesc(squareform(DY)), colormap('gray'), set(gca,'XTick',[1:4]), set(gca,'YTick',[1:4]), axis square, title('ROI2') 


