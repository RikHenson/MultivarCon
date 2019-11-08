clear all
clc
a = figure;
set(a, 'Position', [0, 0, 1200, 350])

nTime = 200; 
%% 1 When Pearson of Mean works
X=randn(1,nTime);
X(2,:)=X(1,:)+0.3*randn(1,nTime);
T=eye(2);
Y=T*X+0.3*randn(2,nTime)+2.5;

abs(corr(mean(X,1)',mean(Y,1)'))

subplot(1,3,1)
scatter(X(1,:),X(2,:),5,'filled')
xlabel('Vox 1')
hold on
scatter(Y(1,:),Y(2,:),5,'filled')
ylabel('Vox 2')
axis([-5 5 -5 5])
set(gca,'FontSize',15)
legend('ROI 1', 'ROI 2','Location','northwest')

%% When Pearson-SVD works
X=randn(1,nTime);
X(2,:)=X(1,:)+0.3*randn(1,nTime);
T=[1,0;0, -1];
Y=T*X+0.3*randn(2,nTime)+2.5;

abs(corr(mean(X,1)',mean(Y,1)'))
abs(corr(mean(X,1)',diff(Y,1)')) %equivalent to SVD, ie projecting along y=-x

subplot(1,3,2)
scatter(X(1,:),X(2,:),5,'filled')
xlabel('Vox 1')
hold on
scatter(Y(1,:),Y(2,:),5,'filled')
ylabel('vox 2')
axis([-5 5 -5 5])
set(gca,'FontSize',15)
legend('ROI 1', 'ROI 2','Location','northwest')

%% When MV is better, because data not dominated by one direction
X=randn(2,nTime);
T=randn(2);
Y=T*X+0.3*randn(2,nTime)+2.5;

% also the UV-metrics seem to work, but it only because of the general
% low dimensionality (dim=2). For high dimensions, only MV-methods work
abs(corr(mean(X,1)',mean(Y,1)'))
abs(corr(diff(X,1)',diff(Y,1)'))
DX = pdist(X','euclidean');DY = pdist(Y','euclidean');corr(DX',DY')

subplot(1,3,3)
scatter(X(1,:),X(2,:),5,'filled')
xlabel('Vox 1')
hold on
scatter(Y(1,:),Y(2,:),5,'filled')
ylabel('vox 2')
axis([-5 5 -5 5])
set(gca,'FontSize',15)
legend('ROI 1', 'ROI 2','Location','northwest')