
clear

SNR = 1;

S.fnum = 11;
S.categories = 0;
S.cat_bound = 0;
S.smoothness = 5;
S.timepoints = 1;

cm = colormap('jet');
N = size(cm,1);


%% All positive
C1 = [1 1/2; 1/2 1];
%C2 = [1 1/2; 1/2 1];
%C12 = C1/2;
%C12 = ones(2)/2;
%C = [C1 C12; C12 C2];
%S.C = C;
%W = [1 1/2; 1/2 1];
%W = [1 1; 1 1];
W = [1 0; 0 1];

rng('default')
y = mvnrnd(zeros(N,2),C1);
y(:,3:4) = y(:,1:2)*W;
y = y + randn(N,4)/SNR;
S.y = y;

S.cc = [0.5 0.5 0 0; 0 0 0.5 0.5];   % average across voxels

[y,cy,fnum] = plot_2reg_2vox(S);

eval(sprintf('print -f%d -dtiff cov_pos.tiff',fnum.Number))
eval(sprintf('print -f%d -dtiff scat_pos.tiff',fnum.Number+1))
eval(sprintf('print -f%d -dtiff time_pos.tiff',fnum.Number+2))

%% No between-region correlation
% C1 = [1 1/2; 1/2 1];
% %C2 = [1 1/2; 1/2 1];
% %C12 = ones(2)*0;
% %S.C = [C1 C12; C12 C2];
% W = zeros(2);
% 
% rng('default')
% y = mvnrnd(zeros(N,2),C1);
% y(:,3:4) = y(:,1:2)*W;
% y = y + randn(N,4)/SNR;
% S.y = y;
% 
% S.cc = [0.5 0.5 0 0; 0 0 0.5 0.5];   % average across voxels
% 
% [y,cy,fnum] = plot_2reg_2vox(S);
% 
% eval(sprintf('print -f%d -dtiff cov_none.tiff',fnum.Number))
% eval(sprintf('print -f%d -dtiff scat_none.tiff',fnum.Number+1))
% eval(sprintf('print -f%d -dtiff time_none.tiff',fnum.Number+2))

%% LDA
%C1 = [1 -1; -1 1];
C1 = [1 1/2; 1/2 1];
%C2 = [1 1; 1 1];
%C2 = [1 1/2; 1/2 1];
%C12 = C1/2;
%C12 = [1 -1; -1 1]/2;
%C12 = 3*[1 0; 0 -1]/4;
%S.C = [C1 C12; C12 C2];
S.categories = 0;
W = [1 0; 0 -1];

rng('default')
y = mvnrnd(zeros(N,2),C1);
y(:,3:4) = y(:,1:2)*W;
y = y + randn(N,4)/SNR;
S.y = y;

S.cc = [0.5 0.5 0 0; 0 0 0.5 0.5];   % average across voxels

[y,cy,fnum] = plot_2reg_2vox(S);

eval(sprintf('print -f%d -dtiff cov_neg.tiff',fnum.Number))
eval(sprintf('print -f%d -dtiff scat_neg.tiff',fnum.Number+1))
eval(sprintf('print -f%d -dtiff time_neg.tiff',fnum.Number+2))

S.categories = 1;
S.cat_bound = 0;
S.cc = [0.5 0.5 0 0; 0 0 0.5 -0.5];  % Projection to LDA

[y,cy,fnum] = plot_2reg_2vox(S);

%eval(sprintf('print -f%d -dtiff cov_neg_lda.tiff',fnum.Number))
eval(sprintf('print -f%d -dtiff scat_neg_lda.tiff',fnum.Number+1))
eval(sprintf('print -f%d -dtiff time_neg_lda.tiff',fnum.Number+2))

%% SVD
S.categories = -1;

ty=[];sy=[];
for r=1:2
    [u,s,v] = svd(y(:,[1:2]+(r-1)*2));
    ty(:,r) = u(:,1);
    sy(:,r) = v(:,1);
end
%[R,p]=corrcoef(ty)

S.cc = -[sy(:,1)' 0 0; 0 0 sy(:,2)'];

[y,cy,fnum] = plot_2reg_2vox(S);

%eval(sprintf('print -f%d -dtiff cov_neg_svd.tiff',fnum.Number))
eval(sprintf('print -f%d -dtiff scat_neg_svd.tiff',fnum.Number+1))
eval(sprintf('print -f%d -dtiff time_neg_svd.tiff',fnum.Number+2))


%% DistCor

C1 = [1 1/2; 1/2 1];
W = [1 0; 0 -1];
SNR = 10;
rng('default')
y = mvnrnd(zeros(N,2),C1);
y(:,3:4) = abs(y(:,1:2));
y = y + randn(N,4)/SNR;
S.y = y;

S.cc = [0.5 0.5 0 0; 0 0 0.5 0.5];

[y,cy,fnum] = plot_2reg_2vox(S);

%R = dcor_dc(S.y(:,1:2),S.y(:,3:4))
[R,T,p,df] = dcor_uc(S.y(:,1:2),S.y(:,3:4))

eval(sprintf('print -f%d -dtiff cov_neg_abs_dc.tiff',fnum.Number))
eval(sprintf('print -f%d -dtiff scat_neg_abs_dc.tiff',fnum.Number+1))
eval(sprintf('print -f%d -dtiff time_neg_abs_dc.tiff',fnum.Number+2))

S.smoothness = 0;
S.categories = 0;
S.SNR = Inf;
%S.cm = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1];
S.y = [
    1 1         1 3
    3 3         3 1 
    4 3         3 2     
    1 2         0 3
	1.1 1.1         1.1 3.1
    3 2         0 1
    ];    
S.cm = repmat([1:size(S.y)]'/(size(S.y,1)+1),1,3);

S.cc = [0.5 0.5 0 0; 0 0 0.5 0.5];

[y,cy,fnum] = plot_2reg_2vox(S);

%R = dcor_dc(S.y(:,1:2),S.y(:,3:4))
[R,T,p,df] = dcor_uc(S.y(:,1:2),S.y(:,3:4))

figure(fnum.Number+1)
subplot(1,2,1), axis([0 4 0 4])
subplot(1,2,2), axis([0 4 0 4])

eval(sprintf('print -f%d -dtiff cov_neg_dc.tiff',fnum.Number))
eval(sprintf('print -f%d -dtiff scat_neg_dc.tiff',fnum.Number+1))
eval(sprintf('print -f%d -dtiff time_neg_dc.tiff',fnum.Number+2))


%% RSA
cm = colormap('jet');
N = size(cm,1)
S.cm = cm;
C1 = [1 1/2; 1/2 1];
S.categories = 1;
S.cat_bound = 0.33;
S.smoothness = 0;
S.timepoints = 0;
SNR = 1;

W = {};
W{1} = [0 0; 0 0];
W{2} = [1 0; 0 -1];
W{3} = [0 0; 0 0];

rng('default')
for t = 1:length(W)
    y = mvnrnd(kron([-10 0 10]',ones(round(N/3),2)),C1);
    %y = mvnrnd(kron([-10 0 10; -10 10 10]',ones(round(N/3),1)),C1);
    y(:,3:4) = y(:,1:2)*W{t};
    y = y + randn(size(y,1),4)/SNR;
    S.y = y;
    
    S.cc = [0.5 0.5 0 0; 0 0 0.5 0.5];   % average across voxels
    
    [y,cy,fnum] = plot_2reg_2vox(S);
    
    rdm = [];
    for r=1:size(y,1);
        for c=(r+1):size(y,1)
            %        rdm(r,c) = corr(y(r,1:2)',y(c,1:2)');
            rdm1(r,c) = sqrt(sum((y(r,1:2) - y(c,1:2)).^2));
            rdm2(r,c) = sqrt(sum((y(r,3:4) - y(c,3:4)).^2));
        end
    end
    
    t1 = triu(rdm1,1); t1 = t1(find(t1>0));
    t2 = triu(rdm2,1); t2 = t2(find(t2>0));
    [R,p] = corr(t1,t2)
    
    eval(sprintf('print -f%d -dtiff cov_neg_rsa%d.tiff',fnum.Number,t))
    eval(sprintf('print -f%d -dtiff scat_neg_rsa%d.tiff',fnum.Number+1,t))
    eval(sprintf('print -f%d -dtiff time_neg_rsa%d.tiff',fnum.Number+2,t))   

    figure(fnum.Number+3),subplot(1,2,1),imagesc(rdm1),axis square, subplot(1,2,2),imagesc(rdm2), axis square
    eval(sprintf('print -f%d -dtiff rdm_rsa%d.tiff',fnum.Number+3,t))   
end

fnum = figure; hold on
line([0 1]',[0 1]')
line([1 2]',[1 0]')
axis([0 2 0 1]); 
set(gca,'YTick',[0 1])
set(gca,'XTick',[0 1 2])

eval(sprintf('print -f%d -dtiff rdm_time.tiff',fnum.Number))


