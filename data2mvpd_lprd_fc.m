function [mvpd,lprd,fc,fc_svd]=data2mvpd_lprd_fc(X,Y,options);
% calculates the MultiVariate Pattern Dependence (MVPD) between two
% multivariate time series (Anzellotti et al. 2017, Plos Comput Biol), the
% linearly predicted representational dissimilarity (LPRD) metric (Basti et al. 2019), 
% the Pearson correlation between the average time series and the one after a singular value decomposition.
%
% Input:
% X and Y:    two cell arrays. The number of cells represents the number of runs,
%             and the dimension of each cell is equal to ntxna, for the ROI1,
%             and ntxnb, for the ROI2.
% options     method: full, svd_exvar or svd_ndir; full denotes taking the whole time
%             series without applying singular value decomposition (SVD). In this case,
%             ridge regression approach is used to estimate the transformation.
%             svd_exvar and svd_ndir are associated with SVD and denote the percentage
%             of variance that we want to explain by applying SVD, or the number
%             of components to consider. In this case, the ordinary least squares
%             approach is used to estimate tha transformation.
% Output:
% mvpd:       MVPD value.
% lprd:       LPRD value (correlation between estimated and actual RDM)
% fc:         Pearson correlation coefficient between average time series.
% fc_svd:     Pearson correlation coefficient between first two SVs.
% Alessio Basti
% version: 29/07/2019

for r=1:length(X)
    % get the voxel-average time-series
    ts_a{r} = mean(X{r},2);
    ts_b{r} = mean(Y{r},2);
    % compute the pearson correlation
    fc_app(r) = corr(ts_a{r},ts_b{r});
    
    opt = options;
    opt.number = 1; % override user arguments for this special case of using SVD to reduce dimension to 1
    opt.meancorrection = 0; % so can return mean over voxels if dominant spatial mode
    [C1_a]=dimreduction(X{r},'svd_ndir',opt);
    [C1_b]=dimreduction(Y{r},'svd_ndir',opt);
    fc_SVs_app(r) = abs(corr(C1_a,C1_b)); % abs because sign of first SV arbitrary
end

method=zeros(2,1);
for jmet=1:1 % do we need this for loop? can't we just do jmet = 1?
    if(jmet==1)
        X_app=X;
        Y_app=Y;
    else
        X_app=ts_a;
        Y_app=ts_b;
    end
    for r=1:length(X_app)
        
        % let us divide the training set from the testing set
        Xtrain=X_app;
        Ytrain=Y_app;
        Xtrain{r}=[];
        Ytrain{r}=[];
        Xtrain=vertcat(Xtrain{:});
        Ytrain=vertcat(Ytrain{:});
        Xtest=X_app{r};
        Ytest=Y_app{r};
        
        if strcmp(options.method,'full')==1
            % linear model estimate (ridge regression)
            [Ttilde,~]=ridgeregmethod(Xtrain,Ytrain,options.regularisation);
            
            % OLS on full time series is not recommended when the number of time
            % points is lower than the number of voxels in X
            %Ttilde=Ytrain'*Xtrain*inv(Xtrain'*Xtrain);
            
            % correlation between the forecasted and the test data
            Ytest_for=Xtest*Ttilde';
            for icomp=1:length(Ytest_for(1,:))
                M=corrcoef(Ytest(:,icomp),Ytest_for(:,icomp));
                method(jmet)=method(jmet)+(var(Ytest(:,icomp))/sum(var(Ytest)))*M(1,2);
            end
        else
            % application of the dimensionality reduction, e.g. selection of the directions which explain
            % a sufficient amount of variance (coded by the input parameter 'percentage')
            [Xtrain_red,Va,SVa]=dimreduction(Xtrain,options.method,options);
            [Ytrain_red,Vb,SVb]=dimreduction(Ytrain,options.method,options);
            
            % linear model estimate (OLS)
            Ttilde=Ytrain_red'*Xtrain_red*inv(Xtrain_red'*Xtrain_red);
            Xtest_red=(Xtest-repmat(mean(Xtest),length(Xtest(:,1)),1))*Va;
            Ytest_red=(Ytest-repmat(mean(Ytest),length(Ytest(:,1)),1))*Vb;
            
            % correlation between the forecasted and the test data
            Ytest_for_red=Xtest_red*Ttilde';
            for icomp=1:length(Ytest_for_red(1,:))
                M=corrcoef(Ytest_red(:,icomp),Ytest_for_red(:,icomp));
                method(jmet)=method(jmet)+(SVb(icomp)/sum(SVb))*M(1,2);
            end
        end
        
        % linear model estimate (ridge regression)
        zX{1}=zscore(X{r},0,2);
        zY{1}=zscore(Y{r},0,2);
        [Ttilde,~]=ridgeregmethod(zX{1},zY{1},options.regularisation);
        zY_for{1}=zX{1}*Ttilde';
        
        % correlation between the estimated and the actual RDM for the ROI2
        [lprd(r,jmet),~] = data2rc(zY,zY_for,'Correlation');
    end
end

mvpd=method(1)/length(X_app);
lprd=mean(lprd(:,1));
fc=mean(fc_app);
fc_svd=mean(fc_SVs_app);

end
