function [mvpd,fc,fc_pc]=data2mvpd(Ya,Yb,options);
% calculates the MultiVariate Pattern Dependence (MVPD) between two
% multivariate time series (Anzellotti et al. 2017, Plos Comput Biol).
%
% Input:
% Ya and Yb:  two cell arrays. The number of cells represents the number of runs, and the dimension of each cell is 
%             equal to ntxna, for the ROIa, and ntxnb, for the ROIb.
% options     method: either pca_exvar or pca_ndir; either percentage or number: the percentage of variance that we want to explain by
%             applying dimensionality reduction approach or the number of
%             ICs to be considered.
% Output: 
% mvpd:       MVPD value.
% uvpd:       UVPD value.
% fc:         Pearson correlation coefficient.
% Alessio Basti 
% version: 04/07/2019

for irun=1:length(Ya)
    %Ya_zs{irun} = zscore(Ya{irun},0,2);
    %Yb_zs{irun} = zscore(Yb{irun},0,2);
    Ya_zs{irun} = Ya{irun};
    Yb_zs{irun} = Yb{irun};
    % get the voxel-average ts
    ts_a{irun} = mean(Ya{irun},2);
    ts_b{irun} = mean(Yb{irun},2);
    % compute the pearson correlation
    fc_app(irun) = corr(ts_a{irun},ts_b{irun});
    [PC1_a{irun}]=dimreduction(Ya_zs{irun},'pca_ndir',options);
    [PC1_b{irun}]=dimreduction(Yb_zs{irun},'pca_ndir',options);
    fc_PCs_app(irun) = corr(PC1_a{irun},PC1_b{irun});
end

method=zeros(2,1);
for jmet=1:2
    if(jmet==1)
       Ya_app=Ya_zs; 
       Yb_app=Yb_zs;
    else
       Ya_app=ts_a; 
       Yb_app=ts_b;
    end
    for irun=1:length(Ya_app)

        %let us divide the training set from the testing set
        Yatrain=Ya_app;
        Ybtrain=Yb_app;
        Yatrain{irun}=[];
        Ybtrain{irun}=[];
        Yatrain=vertcat(Yatrain{:});
        Ybtrain=vertcat(Ybtrain{:});
        Yatest=Ya_app{irun};
        Ybtest=Yb_app{irun};

        % application of the dimensionality reduction, e.g. selection of the directions which explain
        % a sufficient amount of variance (coded by the input parameter 'percentage')
        [Yatrain_red,Va,SVa]=dimreduction(Yatrain,options.method,options);
        [Ybtrain_red,Vb,SVb]=dimreduction(Ybtrain,options.method,options);

        %linear model estimate (least-square method)
        B=Ybtrain_red'*Yatrain_red*pinv(Yatrain_red'*Yatrain_red);
        Yatest_red=(Yatest-repmat(mean(Yatest),length(Yatest(:,1)),1))*Va;
        Ybtest_red=(Ybtest-repmat(mean(Ybtest),length(Ybtest(:,1)),1))*Vb;

        % correlation between the forecasted and the test data
        Ybtest_for_red=Yatest_red*B';
        for icomp=1:length(Ybtest_for_red(1,:))
            M=corrcoef(Ybtest_red(:,icomp),Ybtest_for_red(:,icomp));
            method(jmet)=method(jmet)+(SVb(icomp)/sum(SVb))*M(1,2);
        end
    end
end

%averaging across the runs the MVPD (method(1)) and UVPD (method(2))
mvpd=method(1)/length(Ya_app);
uvpd=method(2)/length(Ya_app);
fc=mean(fc_app);
fc_pc=mean(fc_PCs_app);


end
