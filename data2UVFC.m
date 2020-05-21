function [fc,fc_svd,fc_cca]=data2UVFC(X,Y)
% calculates the Pearson correlation between the average time series and 
% the one after a singular value decomposition.
%
% Input:
% X and Y:    two cell arrays. The number of cells represents the number of runs,
%             and the dimension of each cell is equal to ntxna, for the ROI1,
%             and ntxnb, for the ROI2.
% Output:
% fc:         Pearson correlation coefficient between average time series.
% fc_svd:     Pearson correlation coefficient between first two SVs.
% Hamed
% version: 29/04/2020

for r=1:length(X)
    % get the voxel-average time-series
    ts_a{r} = mean(X{r},2);
    ts_b{r} = mean(Y{r},2);
    % compute the pearson correlation
    fc_app(r) = corr(ts_a{r},ts_b{r});
    
    opt.number = 1; % override user arguments for this special case of using SVD to reduce dimension to 1
    opt.meancorrection = 0; % so can return mean over voxels if dominant spatial mode
    [C1_a]=dimreduction(X{r},'svd_ndir',opt);
    [C1_b]=dimreduction(Y{r},'svd_ndir',opt);
    fc_SVs_app(r) = abs(corr(C1_a,C1_b)); % abs because sign of first SV arbitrary
       
    % Still need to apply some SVD for CCA to work (though still >1 dim hopefully!)
    opt2.percentage = 90; opt2.meancorrection = 1;
    X{r} = dimreduction(X{r},'svd_exvar',opt2);
    Y{r} = dimreduction(Y{r},'svd_exvar',opt2);
              
    % Also check that CCA is estimable (more timepoints than sum of voxels of two ROIs)             
    numTP = size(X{r},1); % Y{r} must have same 
    minTP = size(X{r},2) + size(Y{r},2);
    cca_opt = opt;
    if numTP <= minTP
       cca_opt.number = min((size(Xr{r},2)-1),floor((numTP-1)/2)); % assuming numTP>2!
       Xr{r} = Xr{r}(:,1:cca_opt.number);
       cca_opt.number = min((size(Yr{r},2)-1),floor((numTP-1)/2)); % assuming numTP>2!
       Yr{r} = Yr{r}(:,1:cca_opt.number);
    end    
    [Acca Bcca Rcca]=canoncorr(X{r},Y{r});
    fc_CCA_app(r) = Rcca(1);
end
fc=mean(fc_app);
fc_svd=mean(fc_SVs_app);
fc_cca=mean(fc_CCA_app);

