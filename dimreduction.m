function [data,V,SV]=dimreduction(data,method,options)
% calculates the dimensionality reducted time series by using either the
% singular value decomposition or the averaging process across the voxels.
% input: data:       
%             matrix of dimension Tx(n), where T is the number of
%             time-points or conditions and n is the number of voxels.
% method:     
%             'average' is for averaging across voxels, 'svd_ndir' is
%             for the first N svd directions while 'svd_exvar' is for the
%             directions which allow to explain the percentage of variance
%             indicated by the input parameter (percentage).
% options:    
%             percentage or number. Percentage is in the range [0,100],
%             e.g. 90 is for explaining the 90% of variance. Number is an
%             integer in the range [1,n] (options.percentage needs to be
%             set when method is set to 'svd_ndir' and options.percentage
%             needs to be set for 'svd_exvar') if options.meancorrection is
%             set to 1 then the input data would be centered
%             (mean-corrected across the second dimension, e.g. voxel)
%             before pca is applied.
% Alessio Basti version: 29/07/2019

if ~isfield(options,'meancorrection')
    options.meancorrection  = 1;% overwritten to 0 for PCA FC
end

if strcmp(method,'svd_ndir')==1 
    
    [V newdata SV] = pca(data,'Centered',options.meancorrection);
    
    number=min([options.number; length(data(1,:))]);
    data=newdata(:,1:number);
    V=V(:,1:number);
%    if sign(mean(V)) ~= sign(mean(mean(data))), data = -data; end % reverse sign to match mean
    SV=SV(1:number);
     
elseif strcmp(method,'svd_exvar')==1
     
    [V newdata SV] = pca(data,'Centered',options.meancorrection);
    index = find(cumsum(SV)/sum(SV)>=options.percentage/100,1);
    data=newdata(:,1:index);
    V=V(:,1:index);
    SV=SV(1:index); 
    
elseif strcmp(method,'average')==1  
    
    data=mean(data,2);
    V=[];
    SV=[];
    
end

end
