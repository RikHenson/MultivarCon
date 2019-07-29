function [data,V,SV]=dimreduction(data,method,options);
% calculates the dimensionality reducted time series by using either the
% singular value decomposition or the averaging process across the voxels.
% input:
% data:       matrix of dimension Tx(n)
% method:     'average' is for averaging across voxels, 'svd_ndir' is 
%             for the first N svd directions while 'svd_exvar' is for the directions 
%             which allow to explain the percentage of variance indicated by
%             the input parameter (percentage).
% options:    percentage or number. Percentage is in the range [0,100], e.g. 90 is for explaining the
%             90% of variance. Number is an integer in the range [1,n]
% Alessio Basti 
% version: 29/07/2019

if strcmp(method,'svd_ndir')==1
    
    if(isfield(options,'meancorrection'))
       data=data-repmat(mean(data),length(data(:,1)),1);
    end
    [U S V]=svd(data);
    number=min([options.number; length(data(1,:))]);
    data=data*V(:,1:number);
    V=V(:,1:number);
    SV=var(data);
    
elseif strcmp(method,'svd_exvar')==1
     
    if(isfield(options,'meancorrection'))
       data=data-repmat(mean(data),length(data(:,1)),1);
    end
    [U S V]=svd(data);
    totvar=sum(var(data));
    index=0;
    i=1;
    while index==0
        varapp(i)=var(data*V(:,i));
        sumvar(i)=sum(varapp(1:i));
        if(sumvar(i)>(options.percentage/100)*totvar)
            index=i;
        end
        i=i+1;
    end
    data=data*V(:,1:index);
    V=V(:,1:index);
    SV=var(data);
    
elseif strcmp(method,'average')==1  
    
    data=mean(data,2);
    V=[];
    SV=[];
    
end

end
