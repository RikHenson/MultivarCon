
cd /home/rh01/GitHub/MultivarCon

%% 1 When Pearson of Mean works
d = 0.2;

X = [
    -3  2
    -2  2
    3 -2
    2 -2
    ];

%X = X + d*[0 -1; 0 1; 0 1; 0 -1]; %  + 0.1*randn(4,2);               
Y = X; % + 0.1*randn(4,2);       

corr(mean(X,2),mean(Y,2))

figure, 
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
print -dpng 'UVmean.png'


%% When Pearson of SVD works
X = [
    -3   2
    -1   2
    -3   1
    -1   1
    ];  
Y = [
    1   -2
    3   -2
    1   -1
    3   -1
    ];  
%X = X + d*[-1 1; 1 -1; 1 -1; -1 1]; % + 0.1*randn(4,2);         
%Y = [-X(:,1) X(:,2)];
%Y = Y + d*[-1 -1; 1 1; 1 1; -1 -1]; % + 0.1*randn(4,2);          
corr(mean(X,2),mean(Y,2))
corr(diff(X,1,2),sum(Y,2)) %equivalent to SVD, ie projecting along y=-x

figure, 
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
print -dpng 'UVsvd.png'


%% When MV needed, because data not dominated by one direction (though SVD=diff also works in this case?)

d = 1;

X = [
    2   2
    2   2
    -2  -2
    -2  -2
    ];    
X = X + d*[-1 1; 1 -1; 1 -1; -1 1];% + 0.1*randn(4,2);             
Y = [-X(:,1) X(:,2)];
Y = Y + d*[-1 -1; 1 1; 1 1; -1 -1]; % + 0.1*randn(4,2);               
corr(mean(X,2),mean(Y,2))
corr(diff(X,1,2),sum(Y,2)) % Still -1...

DX = pdist(X,'euclidean');
DY = pdist(Y,'euclidean');
corr(DX',DY')

figure,
subplot(1,2,1); title('ROI1'), imagesc(squareform(DX)), colormap('gray'), axis off
subplot(1,2,2); title('ROI2'), imagesc(squareform(DY)), colormap('gray'), axis off

figure, 
subplot(1,2,1); title('ROI1')
hold on, for t=1:size(X,1); t = text(X(t,1),X(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
subplot(1,2,2); title('ROI2')
hold on, for t=1:size(X,1); t = text(Y(t,1),Y(t,2),num2str(t),'FontSize',12); end; axis square
axis([-4 4 -4 4]), set(gca,'XAxisLocation','origin','YAxisLocation','origin','XTick',[],'YTick',[]), 
print -dpng 'MV1.png'


