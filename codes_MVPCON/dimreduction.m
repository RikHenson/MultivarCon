function [data,V,SV]=dimreduction(data,method,options);
% calculates the dimenionality reducted time series by using either the
% principal component analysis on the covariance matrix or the averaging
% process across the voxels.
% input:
% data:       matrix of dimension Tx(n)
% method:     'average' is for averaging across voxels, 'pca_dir' is 
%             for the first N PCA directions while 'pca' is for the directions 
%             which allow to explain the percentage of variance indicated by ù
%             the input parameter (percentage).
% options:    percentage or number. Percentage is in the range [0,100], e.g. 90 is for explaining the
%             90% of variance. Number is an integer in the range [1,n]
% Alessio Basti 
% version: 12/10/2018

if strcmp(method,'pca_ndir')==1
    
    [V newdata SV]=pca(data-repmat(mean(data),length(data(:,1)),1));
    number=min([options.number; length(data(1,:))]);
    data=newdata(:,1:number);
    V=V(:,1:number);
    SV=SV(1:number);
    
elseif strcmp(method,'pca_exvar')==1
       
    [V newdata SV]=pca(data-repmat(mean(data),length(data(:,1)),1));
    index=0;
    i=1;
    while index==0
        varapp(i)=var(newdata(:,i));
        sumvar(i)=sum(varapp(1:i));
        if(sumvar(i)>(options.percentage/100)*sum(SV))
            index=i;
        end
        i=i+1;
    end
    data=newdata(:,1:index);
    V=V(:,1:index);
    SV=SV(1:index); 
    
elseif strcmp(method,'average')==1  
    
    data=mean(data,2);
    V=[];
    SV=[];
    
end

end